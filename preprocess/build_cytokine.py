"""Build the cytokine dictionary stores (mouse Immune Dictionary, human Cytokine Dictionary).

    python -m preprocess.build_cytokine --species mouse|human

Writes app_data/cytokine_<species>/:
    logFC.npy                        float32 [n_tested, n_celltypes, n_cytokines]   NaN = not tested in that cell type
    PValue.npy, FDR.npy              float64, same shape (float32 underflows below ~1e-45)
    pct_expr.npy                     float32 [n_genes, n_celltypes, n_cytokines + 1] last slot = PBS
    n_cells.npy                      int64   [n_celltypes, n_cytokines + 1]
    genes.parquet, celltypes.parquet, cytokines.parquet, meta.json
"""
from __future__ import annotations

import argparse
import io
import json
import re
import tarfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from preprocess.annotations import read_gtf_genes, read_h5ad_obs_column, read_h5ad_var
from scbrowser.config import app_data, get_path, load_config, resolve
from scbrowser.ids import strip_ens

METRICS = ("logFC", "PValue", "FDR")
DTYPES = {"logFC": np.float32, "PValue": np.float64, "FDR": np.float64}

# Order and display labels follow the manuscript dot plot (Figure_examples/dotplotexample.png).
MOUSE_CELLTYPES = [
    ("B_cell", "B cell"), ("T_cell_CD4", "CD4 T"), ("Treg", "Treg"), ("T_cell_CD8", "CD8 T"),
    ("T_cell_gd", "γδ T"), ("NK_cell", "NK"), ("cDC1", "cDC1"), ("cDC2", "cDC2"), ("MigDC", "MigDC"),
    ("Langerhans", "Langerhans"), ("pDC", "pDC"), ("Macrophage", "Macrophage"),
]
HUMAN_UNMAPPED_CELLTYPES = ["HSPC", "ILC", "Plasmablast"]   # no mouse counterpart, not in the lineage map

FAMILY_COLUMNS = ["DEG_plot_family", "single_family", "family_level_1", "family_level_2"]
FAMILY_SHORT = {
    "Common γ chain": "γc", "Common β chain": "βc", "Interferon": "IFN", "Growth factor": "Growth",
    "Complement": "Compl.", "IL-4/IL-13": "IL-4/13", "interferon": "IFN", "growth factor": "Growth",
    "complement": "Compl.", "interleukin": "IL", "TNF family": "TNF", "chemokine": "Chemokine",
}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def read_tsv_member(tar: tarfile.TarFile, name: str, **kw) -> pd.DataFrame:
    raw = tar.extractfile(name).read()
    return pd.read_csv(io.BytesIO(raw), sep="\t", compression="gzip" if name.endswith(".gz") else None, **kw)


def file_label(display: str) -> str:
    return display.replace(" ", "_").replace("/", "_")


def celltype_table(species: str, deposit: Path) -> pd.DataFrame:
    lineage = pd.read_csv(deposit / "celltype_lineage_map.tsv", sep="\t")
    lineage = lineage[lineage.species == species]
    if species == "mouse":
        lin = dict(zip(lineage.cell_type, lineage.lineage))
        rows = [(label, display, lin.get(label, "")) for label, display in MOUSE_CELLTYPES]
    else:
        rows = [(file_label(ct), ct, lin) for ct, lin in zip(lineage.cell_type, lineage.lineage)]
        rows += [(ct, ct, "") for ct in HUMAN_UNMAPPED_CELLTYPES]
    return pd.DataFrame(rows, columns=["label", "display", "lineage"])


def cytokine_table(species: str, labels: list[str], deposit: Path) -> pd.DataFrame:
    """One row per tested cytokine, in the authors' order, with display names and family columns."""
    fam = pd.read_csv(resolve(load_config()["cytokine_families"]), sep="\t").set_index("h5ad_label")
    ligand_map = {}
    if species == "mouse":
        lm = pd.read_csv(deposit / "ligand_name_map.tsv", sep="\t")
        ligand_map = dict(zip(lm.mouse_ligand, lm.human_cytokine))
    rows, unmapped = [], []
    for label in labels:
        human_label = ligand_map.get(label, label)
        if human_label in fam.index:
            r = fam.loc[human_label]
            display = label if species == "mouse" else r["authors_name"]
            rows.append(dict(label=label, display=display, order=int(r["order"]), mapped=True,
                             **{c: r[c] for c in FAMILY_COLUMNS}))
        else:
            unmapped.append(label)
    # cytokines the family table doesn't know (mouse-only) go to "Other", after the listed ones
    for k, label in enumerate(sorted(unmapped, key=str.lower)):
        rows.append(dict(label=label, display=label, order=1000 + k, mapped=False,
                         DEG_plot_family="Other", single_family="other", family_level_1="other",
                         family_level_2="other"))
    if unmapped:
        log(f"cytokines without a family (assigned to Other): {', '.join(unmapped)}")
    df = pd.DataFrame(rows).sort_values("order", kind="stable").reset_index(drop=True)
    for c in FAMILY_COLUMNS:
        df[c + "_short"] = df[c].map(lambda v: FAMILY_SHORT.get(v, v))
    return df


def fill_de(arrays: dict[str, np.ndarray], c: int, metric: str, table: pd.DataFrame,
            row_of: pd.Series, col_of: dict[str, int]) -> None:
    rows = row_of.loc[table.index].to_numpy()
    cols = np.array([col_of[x] for x in table.columns])
    arrays[metric][np.ix_(rows, [c], cols)] = table.to_numpy(arrays[metric].dtype)[:, None, :]


# -- mouse -------------------------------------------------------------------------------------------------------
def build_mouse(out: Path) -> dict:
    deposit = get_path("deposits", "immdict")
    h5ad = deposit / "immdict_mouse_cellbender_counts.h5ad"
    var = read_h5ad_var(h5ad)
    genes = pd.DataFrame({"symbol": var.gene_name, "gene_ids": var.gene_id, "gene_id_version": var.gene_id_version})
    cts = celltype_table("mouse", deposit)

    log("reading mouse genome-wide DE tables")
    tables: dict[tuple[str, str], pd.DataFrame] = {}
    with tarfile.open(deposit / "mouse_edgeR.tar") as tar:
        for ct in cts.label:
            for metric in METRICS:
                df = read_tsv_member(tar, f"mouse_edgeR/genomewide/gw_{ct}_{metric}.tsv.gz")
                tables[ct, metric] = df.drop(columns=["gene_name"]).set_index("gene_id")
        panel = pd.concat([read_tsv_member(tar, f"mouse_edgeR/panel/de_{ct}.tsv") for ct in cts.label])

    labels = sorted({c for (ct, m), df in tables.items() for c in df.columns})
    cyts = cytokine_table("mouse", labels, deposit)
    col_of = {lab: j for j, lab in enumerate(cyts.label)}

    var_row = pd.Series(np.arange(len(genes)), index=genes.gene_id_version)
    tested = sorted({var_row[g] for df in (tables[ct, "logFC"] for ct in cts.label) for g in df.index})
    de_row = np.full(len(genes), -1, np.int32)
    de_row[tested] = np.arange(len(tested))
    row_of = pd.Series(de_row[var_row.to_numpy()], index=var_row.index)

    shape = (len(tested), len(cts), len(cyts))
    arrays = {m: np.full(shape, np.nan, DTYPES[m]) for m in METRICS}
    for c, ct in enumerate(cts.label):
        for metric in METRICS:
            fill_de(arrays, c, metric, tables[ct, metric], row_of, col_of)
    genes["de_row"] = de_row

    log("computing % cells expressing per cell type x condition from the h5ad")
    pct, n_cells = mouse_pct_expressing(h5ad, cts, cyts, len(genes))

    write_store(out, genes, cts, cyts, arrays, pct, n_cells)
    check_mouse(out, genes, cts, cyts, arrays, pct, panel, deposit)
    return {"species": "mouse", "n_genes": len(genes), "n_tested": len(tested)}


def mouse_pct_expressing(h5ad: Path, cts: pd.DataFrame, cyts: pd.DataFrame, n_genes: int):
    ct_cats, ct_codes = read_h5ad_obs_column(h5ad, "cell_type")
    cy_cats, cy_codes = read_h5ad_obs_column(h5ad, "cytokine")
    n_cyt = len(cyts)
    ct_pos = {lab: i for i, lab in enumerate(cts.label)}
    cy_pos = {lab: j for j, lab in enumerate(cyts.label)} | {"PBS": n_cyt}
    row_ct = np.array([ct_pos.get(c, -1) for c in ct_cats])[ct_codes]
    row_cy = np.array([cy_pos.get(c, -1) for c in cy_cats])[cy_codes]
    keep = (row_ct >= 0) & (row_cy >= 0)
    group = np.where(keep, row_ct * (n_cyt + 1) + row_cy, -1)
    n_groups = len(cts) * (n_cyt + 1)
    acc = np.zeros((n_groups, n_genes), np.int64)
    with h5py.File(h5ad, "r") as h:
        X = h["X"]
        indptr = X["indptr"][:].astype(np.int64)
        n_cells_total = len(indptr) - 1
        step = 40_000
        for s in range(0, n_cells_total, step):
            e = min(n_cells_total, s + step)
            a, b = indptr[s], indptr[e]
            # int32 throughout: a bool @ bool sparse product would saturate at 1
            nonzero = (X["data"][a:b] > 0).astype(np.int32)
            block = sp.csr_matrix((nonzero, X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, n_genes))
            g = group[s:e]
            sel = np.flatnonzero(g >= 0)
            onehot = sp.csr_matrix((np.ones(len(sel), np.int32), (g[sel], sel)), shape=(n_groups, e - s))
            acc += (onehot @ block).toarray()
            log(f"  rows {e:,}/{n_cells_total:,}")
    n_cells = np.bincount(group[keep], minlength=n_groups).astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = (acc / n_cells[:, None] * 100).astype(np.float32)
    pct = pct.reshape(len(cts), n_cyt + 1, n_genes).transpose(2, 0, 1)
    return np.ascontiguousarray(pct), n_cells.reshape(len(cts), n_cyt + 1)


def check_mouse(out, genes, cts, cyts, arrays, pct, panel, deposit) -> None:
    """The panel tables must equal the genome-wide ones (README); % expressing must be close to the figure tables."""
    row = pd.Series(genes.de_row.to_numpy(), index=genes.gene_id_version)
    c_of = {lab: i for i, lab in enumerate(cts.label)}
    j_of = {lab: j for j, lab in enumerate(cyts.label)}
    r = row.loc[panel.gene_id].to_numpy()
    assert (r >= 0).all(), "panel gene missing from the genome-wide tables"
    c = panel.cell_type.map(c_of).to_numpy()
    j = panel.cytokine.map(j_of).to_numpy()
    d_lfc = np.abs(arrays["logFC"][r, c, j] - panel.logFC.to_numpy())
    rel_p = np.abs(arrays["PValue"][r, c, j] - panel.PValue.to_numpy()) / np.maximum(panel.PValue.to_numpy(), 1e-300)
    log(f"check panel: {len(panel)} values, max |dlogFC| {d_lfc.max():.2g}, max rel dP {rel_p.max():.2g}")
    assert d_lfc.max() < 1e-3 and rel_p.max() < 1e-2

    with tarfile.open(deposit / "figure_source_tables.tar") as tar:
        fig = read_tsv_member(tar, "figure_source_tables/mouse_pct_cells_by_condition.tsv")
    wanted = set(fig.gene)
    sym_row = {s: i for i, s in enumerate(genes.symbol) if s in wanted}
    fig = fig[fig.gene.isin(sym_row) & fig.cell_type.isin(c_of)]
    j_cond = fig.cytokine.map(lambda x: len(cyts) if x == "PBS" else j_of.get(x, -1)).to_numpy()
    fig = fig[j_cond >= 0]
    j_cond = j_cond[j_cond >= 0]
    mine = pct[fig.gene.map(sym_row).to_numpy(), fig.cell_type.map(c_of).to_numpy(), j_cond]
    diff = np.abs(mine - fig.pct_cells.to_numpy())
    diff = diff[~np.isnan(diff)]
    log(f"check % expressing vs figure table: {len(diff)} values, median {np.median(diff):.3f}, "
        f"95% {np.percentile(diff, 95):.2f}, max {diff.max():.2f} points")
    assert diff.max() <= 7.5, "README states recomputed values differ by at most 7.3 points"


# -- human -------------------------------------------------------------------------------------------------------
def build_human(out: Path) -> dict:
    deposit = get_path("deposits", "immdict")
    h5 = deposit / "human_pseudobulk_all.h5"
    with h5py.File(h5, "r") as h:
        symbols = list(h["genes"].asstr()[:])
    genes = human_gene_table(symbols)
    cts = celltype_table("human", deposit)

    log("reading human genome-wide DE tables")
    tables: dict[tuple[str, str], pd.DataFrame] = {}
    with tarfile.open(deposit / "human_edgeR.tar") as tar:
        for ct in cts.label:
            for metric in METRICS:
                df = read_tsv_member(tar, f"human_edgeR/genomewide/gw_{ct}_{metric}.tsv.gz")
                tables[ct, metric] = df.set_index("gene")
        panel = pd.concat([read_tsv_member(tar, f"human_edgeR/panel/erv_de_{ct}.tsv") for ct in cts.label])

    labels = sorted({c for (ct, m), df in tables.items() for c in df.columns})
    cyts = cytokine_table("human", labels, get_path("deposits", "immdict"))
    col_of = {lab: j for j, lab in enumerate(cyts.label)}

    sym_row = pd.Series(np.arange(len(genes)), index=genes.symbol)
    tested = sorted({sym_row[g] for df in (tables[ct, "logFC"] for ct in cts.label) for g in df.index})
    de_row = np.full(len(genes), -1, np.int32)
    de_row[tested] = np.arange(len(tested))
    row_of = pd.Series(de_row[sym_row.to_numpy()], index=sym_row.index)
    shape = (len(tested), len(cts), len(cyts))
    arrays = {m: np.full(shape, np.nan, DTYPES[m]) for m in METRICS}
    for c, ct in enumerate(cts.label):
        for metric in METRICS:
            fill_de(arrays, c, metric, tables[ct, metric], row_of, col_of)
    genes["de_row"] = de_row

    log("summing n_expressing over donors per cell type x condition")
    pct, n_cells, total_expr = human_pct_expressing(h5, cts, cyts, len(genes))

    write_store(out, genes, cts, cyts, arrays, pct, n_cells)
    check_human(genes, cts, cyts, arrays, total_expr, panel, deposit)
    return {"species": "human", "n_genes": len(genes), "n_tested": len(tested)}


def human_gene_table(symbols: list[str]) -> pd.DataFrame:
    v43 = read_gtf_genes(resolve(load_config()["annotations"]["human_v43_gtf"]))
    v43 = v43[~v43.par_y]
    name_ids = v43.groupby("gene_name")["gene_id"].agg(lambda s: sorted(set(s))).to_dict()
    present = set(symbols)
    ids, variant_of = [], []
    for s in symbols:
        base = None
        if s.startswith("ENSG"):
            gids = [strip_ens(s)]
        elif s in name_ids:
            gids = name_ids[s]
        else:
            # anndata make_unique variants (MATR3-1) of a symbol present in the data; real symbols like H1-1 are in v43
            m = re.match(r"^(.+)-(\d+)$", s)
            base = m.group(1) if m and m.group(1) in present else None
            gids = name_ids.get(base, []) if base else []
        ids.append("|".join(gids))
        variant_of.append(base)
    n_var = sum(v is not None for v in variant_of)
    n_none = sum(not i for i in ids)
    log(f"human genes: {len(symbols):,}; make_unique variants {n_var}; without an Ensembl ID {n_none}")
    return pd.DataFrame({"symbol": symbols, "gene_ids": ids, "variant_of": variant_of})


def human_pct_expressing(h5: Path, cts: pd.DataFrame, cyts: pd.DataFrame, n_genes: int):
    n_cyt = len(cyts)
    ct_pos = {d: i for i, d in enumerate(cts.display)}
    cy_pos = {lab: j for j, lab in enumerate(cyts.label)} | {"PBS": n_cyt}
    n_groups = len(cts) * (n_cyt + 1)
    with h5py.File(h5, "r") as h:
        g_ct = h["groups/cell_type"].asstr()[:]
        g_cy = h["groups/cytokine"].asstr()[:]
        g_n = h["groups/n_cells"][:].astype(np.int64)
        missing = (set(g_ct) - set(ct_pos)) | (set(g_cy) - set(cy_pos))
        assert not missing, f"groups with unknown labels: {sorted(missing)[:10]}"
        group = np.array([ct_pos[a] * (n_cyt + 1) + cy_pos[b] for a, b in zip(g_ct, g_cy)])
        acc = np.zeros((n_groups, n_genes), np.int64)
        ne = h["n_expressing"]
        step = 2048
        for s in range(0, ne.shape[0], step):
            e = min(ne.shape[0], s + step)
            onehot = sp.csr_matrix((np.ones(e - s, np.int64), (group[s:e], np.arange(e - s))), shape=(n_groups, e - s))
            acc += onehot @ ne[s:e].astype(np.int64)
            log(f"  groups {e:,}/{ne.shape[0]:,}")
    n_cells = np.bincount(group, weights=g_n, minlength=n_groups).astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = (acc / n_cells[:, None] * 100).astype(np.float32)
    pct = np.ascontiguousarray(pct.reshape(len(cts), n_cyt + 1, n_genes).transpose(2, 0, 1))
    return pct, n_cells.reshape(len(cts), n_cyt + 1), acc.sum(axis=0)


def check_human(genes, cts, cyts, arrays, total_expr, panel, deposit) -> None:
    with tarfile.open(deposit / "figure_source_tables.tar") as tar:
        fig = read_tsv_member(tar, "figure_source_tables/human_pct_cells.tsv")
    mine = pd.Series(total_expr, index=genes.symbol).loc[fig.gene].to_numpy()
    n_bad = int((mine != fig.n_cells_expressing.to_numpy()).sum())
    log(f"check overall expressing-cell counts vs human_pct_cells.tsv: {n_bad} mismatches of {len(fig):,}")
    assert n_bad == 0

    row = pd.Series(genes.de_row.to_numpy(), index=genes.symbol)
    c_of = {lab: i for i, lab in enumerate(cts.label)}
    j_of = {lab: j for j, lab in enumerate(cyts.label)}
    r = row.loc[panel.gene].to_numpy()
    ok = r >= 0
    c = panel.cell_type.map(c_of).to_numpy()[ok]
    j = panel.cytokine.map(j_of).to_numpy()[ok]
    d_lfc = np.abs(arrays["logFC"][r[ok], c, j] - panel.logFC.to_numpy()[ok])
    pv = panel.PValue.to_numpy()[ok]
    rel_p = np.abs(arrays["PValue"][r[ok], c, j] - pv) / np.maximum(pv, 1e-300)
    log(f"check ERV panel: {ok.sum()} values, max |dlogFC| {np.nanmax(d_lfc):.2g}, max rel dP {np.nanmax(rel_p):.2g}")
    assert np.nanmax(d_lfc) < 1e-3 and np.nanmax(rel_p) < 1e-2


# -- output ------------------------------------------------------------------------------------------------------
def write_store(out: Path, genes, cts, cyts, arrays, pct, n_cells) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for m in METRICS:
        np.save(out / f"{m}.npy", arrays[m])
    np.save(out / "pct_expr.npy", pct)
    np.save(out / "n_cells.npy", n_cells)
    genes.to_parquet(out / "genes.parquet", index=False)
    cts.to_parquet(out / "celltypes.parquet", index=False)
    cyts.to_parquet(out / "cytokines.parquet", index=False)
    log(f"wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", choices=["mouse", "human"], required=True)
    args = ap.parse_args()
    out = app_data() / f"cytokine_{args.species}"
    meta = build_mouse(out) if args.species == "mouse" else build_human(out)
    meta["fdr"] = "Benjamini-Hochberg over all genes within each cell type x cytokine (deposit genome-wide tables)"
    (out / "meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
