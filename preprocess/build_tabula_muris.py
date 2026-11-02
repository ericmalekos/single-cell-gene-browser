"""Build the Tabula Muris store from the SCE .rds files.

    python -m preprocess.build_tabula_muris [--keep-export] [--immune-only]

Runs export_tm_sce.R on the global all-ages object (logcounts + colData + UMAP) and, when configured, on the
immune subset (colData + UMAP + a few check genes). Writes app_data/tm/:
    expr/indptr.npy, indices.npy, data.npy   logcounts as gene-major CSR (one gene = one contiguous slice)
    genes.parquet                           symbol, gene_ids (vM38), gene_type, n_features, n_expr, is_spike
    cells_<embedding>.parquet               cell_id, x, y, metadata categoricals[, expr_index]
    embeddings.json                         per embedding: label, cells file, expression dir, reduced dim

The immune subset reuses the global expression store through `expr_index` (its cells joined to global cells by
cell id), after checking that the check genes have identical values. If the join or the check fails, the
subset gets its own store aligned to the global gene order. --immune-only builds just the subset (development
use while the global object is unavailable).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from preprocess.annotations import read_gtf_genes
from scbrowser.config import app_data, get_path, load_config, resolve
from scbrowser.ids import strip_ens

# colData column -> store column
CELL_COLUMNS = {
    "tissue": "tissue", "cell_ontology_class": "cell_type", "free_annotation": "free_annotation", "age": "age",
    "sex": "sex", "mouse.id": "mouse_id", "subtissue": "subtissue", "cluster_louvain": "cluster",
}
EMBEDDING_LABELS = {"global": "All cells, all ages", "immune": "Immune cells, all ages"}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def md5sum(path: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def run_export(rds: Path, outdir: Path, matrix: bool, genes: list[str] = ()) -> None:
    if (outdir / "done").exists():
        log(f"reusing export in {outdir}")
        return
    script = Path(__file__).with_name("export_tm_sce.R")
    cmd = ["Rscript", str(script), "--rds", str(rds), "--outdir", str(outdir)]
    if not matrix:
        cmd.append("--no-matrix")
    if genes:
        cmd += ["--genes", ",".join(genes)]
    log("running " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    (outdir / "done").touch()


def read_cells(d: Path) -> pd.DataFrame:
    raw = pd.read_csv(d / "cells.tsv", sep="\t", dtype=str, keep_default_na=False)
    umap = pd.read_csv(d / "umap.tsv", sep="\t", dtype={"cell_id": str})
    assert (raw.cell_id.to_numpy() == umap.cell_id.to_numpy()).all(), "cells.tsv and umap.tsv disagree"
    cells = pd.DataFrame({"cell_id": raw.cell_id,
                          "x": umap.umap_1.to_numpy(np.float32), "y": umap.umap_2.to_numpy(np.float32)})
    for src, dst in CELL_COLUMNS.items():
        if src not in raw:
            continue
        values = raw[src].replace({"nan": "", "NA": ""})
        if dst == "age":
            values = values.replace({"": "unknown"})   # 661 global cells carry no age
            known = sorted(set(values) - {"unknown"}, key=lambda a: float(a.rstrip("m")))
            levels = known + (["unknown"] if "unknown" in set(values) else [])
            cells[dst] = pd.Categorical(values, categories=levels, ordered=True)
        elif dst == "cluster":
            levels = sorted(set(values), key=lambda c: (len(c), c))
            cells[dst] = pd.Categorical(values, categories=levels)
        else:
            cells[dst] = pd.Categorical(values)
    return cells


def build_expression(d: Path, out: Path, gene_order: list[str] | None = None) -> list[str]:
    """Export slots -> gene-major CSR memmaps. With gene_order, rows are re-indexed to that order."""
    n_genes, n_cells, nnz = map(int, (d / "dims.txt").read_text().split())
    symbols = (d / "genes.txt").read_text().splitlines()
    log(f"loading logcounts: {n_genes:,} genes x {n_cells:,} cells, {nnz:,} non-zeros")
    x = np.fromfile(d / "x.f32", dtype="<f4")
    i = np.fromfile(d / "i.i32", dtype="<i4")
    p = np.fromfile(d / "p.i32", dtype="<i4")
    assert len(x) == len(i) == nnz and len(p) == n_cells + 1 and p[-1] == nnz
    csr = sp.csc_matrix((x, i, p), shape=(n_genes, n_cells)).tocsr()
    del x, i, p
    csr.sort_indices()
    if gene_order is not None:
        pos = {s: k for k, s in enumerate(symbols)}
        dropped = set(symbols) - set(gene_order)
        if dropped:
            log(f"warning: {len(dropped)} genes absent from the target gene order are dropped")
        take = np.array([pos.get(s, -1) for s in gene_order])
        empty = sp.csr_matrix((1, n_cells), dtype=np.float32)
        csr = sp.vstack([empty, csr], format="csr")[take + 1]   # row 0 = empty row for missing genes
        symbols = list(gene_order)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "indptr.npy", csr.indptr.astype(np.int64))
    np.save(out / "indices.npy", csr.indices.astype(np.int32))
    np.save(out / "data.npy", csr.data.astype(np.float32))
    (out / "shape.json").write_text(json.dumps({"n_genes": csr.shape[0], "n_cells": csr.shape[1]}))
    log(f"wrote {out} ({csr.nnz:,} non-zeros)")
    return symbols


def gene_table(symbols: list[str], expr_dir: Path) -> pd.DataFrame:
    m38 = read_gtf_genes(resolve(load_config()["annotations"]["mouse_m38_gtf"]))
    m38 = m38[~m38.par_y]
    by_name = m38.groupby("gene_name").agg(ids=("gene_id", lambda s: sorted(set(s))),
                                           types=("gene_type", lambda s: sorted(set(s))))
    ids, types = [], []
    for s in symbols:
        if s in by_name.index:
            ids.append(by_name.at[s, "ids"])
            types.append("|".join(by_name.at[s, "types"]))
        elif s.upper().startswith("ENSMUSG"):
            ids.append([strip_ens(s)])
            types.append("")
        else:
            ids.append([])
            types.append("spike-in" if s.startswith("ERCC-") else "")
    n_expr = np.diff(np.load(expr_dir / "indptr.npy"))
    genes = pd.DataFrame({
        "symbol": symbols,
        "gene_ids": ["|".join(x) for x in ids],
        "gene_type": types,
        "n_features": [len(x) for x in ids],   # >1: tx2gene summed several GENCODE genes sharing the symbol
        "n_expr": n_expr.astype(np.int64),
        "is_spike": [s.startswith("ERCC-") for s in symbols],
    })
    log(f"genes: {len(genes):,}; with vM38 IDs {sum(bool(x) for x in ids):,}; summed symbols "
        f"{int((genes.n_features > 1).sum())}")
    return genes


def gene_vector(expr_dir: Path, row: int, n_cells: int) -> np.ndarray:
    indptr = np.load(expr_dir / "indptr.npy", mmap_mode="r")
    indices = np.load(expr_dir / "indices.npy", mmap_mode="r")
    data = np.load(expr_dir / "data.npy", mmap_mode="r")
    v = np.zeros(n_cells, np.float32)
    s, e = indptr[row], indptr[row + 1]
    v[indices[s:e]] = data[s:e]
    return v


def immune_matches_global(imm_dir: Path, imm: pd.DataFrame, glob: pd.DataFrame, symbols: list[str],
                          expr_dir: Path) -> bool:
    gpos = pd.Series(np.arange(len(glob)), index=glob.cell_id)
    missing = int((~imm.cell_id.isin(gpos.index)).sum())
    if missing:
        log(f"immune join: {missing:,} of {len(imm):,} cells are not in the global object")
        return False
    imm["expr_index"] = gpos.loc[imm.cell_id].to_numpy(np.int64)
    check = pd.read_csv(imm_dir / "check_genes.tsv", sep="\t", dtype={"cell_id": str})
    row_of = {s: k for k, s in enumerate(symbols)}
    for g in check.columns[1:]:
        if g not in row_of:
            log(f"immune check: {g} not in the global genes")
            return False
        v = gene_vector(expr_dir, row_of[g], len(glob))[imm.expr_index.to_numpy()]
        diff = float(np.abs(v - check[g].to_numpy(np.float32)).max())
        log(f"immune check {g}: max |difference| {diff:.2g}")
        if diff > 1e-4:
            return False
    return True


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--keep-export", action="store_true", help="keep the R export (several GB)")
    ap.add_argument("--immune-only", action="store_true", help="build only the immune subset")
    ap.add_argument("--skip-md5", action="store_true")
    args = ap.parse_args(argv)

    cfg = load_config()["tabula_muris"]
    out = app_data() / "tm"
    export = app_data() / "_tm_export"
    global_rds = get_path("tabula_muris", "global_rds")
    immune_rds = get_path("tabula_muris", "immune_rds")
    if immune_rds is not None and not immune_rds.exists():
        log(f"warning: immune subset {immune_rds} not found; building the global embedding only")
        immune_rds = None
    if args.immune_only and immune_rds is None:
        raise SystemExit("--immune-only needs tabula_muris.immune_rds in the config")
    check_genes = cfg.get("check_genes", [])
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    embeddings = []

    if args.immune_only:
        run_export(immune_rds, export / "immune", matrix=True)
        symbols = build_expression(export / "immune", out / "expr")
        imm = read_cells(export / "immune")
        imm.to_parquet(out / "cells_immune.parquet", index=False)
        embeddings.append(dict(key="immune", label=EMBEDDING_LABELS["immune"], cells="cells_immune.parquet",
                               expr="expr", reduced_dim="UMAP"))
    else:
        if not args.skip_md5:
            log("checking md5 of the global object")
            got = md5sum(global_rds)
            assert got == cfg["global_rds_md5"], f"md5 mismatch for {global_rds}: {got}"
        run_export(global_rds, export / "global", matrix=True)
        symbols = build_expression(export / "global", out / "expr")
        glob = read_cells(export / "global")
        glob.to_parquet(out / "cells_global.parquet", index=False)
        embeddings.append(dict(key="global", label=EMBEDDING_LABELS["global"], cells="cells_global.parquet",
                               expr="expr", reduced_dim="UMAP_refined"))

        if immune_rds is not None:
            run_export(immune_rds, export / "immune", matrix=False, genes=check_genes)
            imm = read_cells(export / "immune")
            if immune_matches_global(export / "immune", imm, glob, symbols, out / "expr"):
                expr_dir = "expr"
            else:
                log("building a separate expression store for the immune subset")
                run_export(immune_rds, export / "immune_matrix", matrix=True)
                build_expression(export / "immune_matrix", out / "expr_immune", gene_order=symbols)
                imm = imm.drop(columns=["expr_index"], errors="ignore")
                expr_dir = "expr_immune"
            imm.to_parquet(out / "cells_immune.parquet", index=False)
            embeddings.append(dict(key="immune", label=EMBEDDING_LABELS["immune"], cells="cells_immune.parquet",
                                   expr=expr_dir, reduced_dim="UMAP"))

    gene_table(symbols, out / "expr").to_parquet(out / "genes.parquet", index=False)
    (out / "embeddings.json").write_text(json.dumps(embeddings, indent=2))
    if not args.keep_export:
        shutil.rmtree(export)
        log("removed the R export")
    log("done")


if __name__ == "__main__":
    main()
