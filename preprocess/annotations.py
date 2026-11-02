"""Gene annotation tables: GTF gene records, h5ad var, symbol aliases (MGI, HGNC) and the cross-dataset ID table.

Run as `python -m preprocess.annotations` to write app_data/annotations/{aliases_auto.tsv, id2symbol.parquet}.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

import h5py
import pandas as pd

from scbrowser.config import app_data, get_path, load_config, resolve
from scbrowser.ids import strip_ens

_ATTR = re.compile(r'(gene_id|gene_name|gene_type) "([^"]*)"')


def ann_dir() -> Path:
    d = app_data() / "annotations"
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_gtf_genes(path: Path) -> pd.DataFrame:
    """Gene records of a GENCODE GTF: gene_id (unversioned), gene_id_version, gene_name, gene_type, par_y.

    Parsed once and cached next to the other annotation tables.
    """
    cache = ann_dir() / (Path(path).name.replace(".gtf.gz", "").replace(".gtf", "") + ".genes.parquet")
    if cache.exists() and cache.stat().st_mtime > Path(path).stat().st_mtime:
        return pd.read_parquet(cache)
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = dict(_ATTR.findall(parts[8]))
            gid = attrs["gene_id"]
            rows.append((strip_ens(gid), gid, attrs.get("gene_name", gid), attrs.get("gene_type", ""),
                         gid.endswith("_PAR_Y")))
    df = pd.DataFrame(rows, columns=["gene_id", "gene_id_version", "gene_name", "gene_type", "par_y"])
    df.to_parquet(cache, index=False)
    return df


def _h5_strings(node) -> list[str]:
    """Read an anndata string column (plain or categorical) as UTF-8 strings."""
    if isinstance(node, h5py.Group):
        cats = node["categories"].asstr()[:]
        return [cats[c] if c >= 0 else "" for c in node["codes"][:]]
    return list(node.asstr()[:])


def read_h5ad_var(path: Path) -> pd.DataFrame:
    """var of an h5ad (no X): gene_id_version (index), gene_id, gene_name, gene_type if present."""
    with h5py.File(path, "r") as h:
        var = h["var"]
        index_key = var.attrs.get("_index", "_index")
        out = pd.DataFrame({"gene_id_version": _h5_strings(var[index_key])})
        for col in ("gene_name", "gene_type", "feature_type"):
            if col in var:
                out[col] = _h5_strings(var[col])
    out["gene_id"] = out["gene_id_version"].map(strip_ens)
    return out


def read_h5ad_obs_column(path: Path, column: str) -> tuple[list[str], "pd.Series"]:
    """Categories and integer codes of a categorical obs column."""
    with h5py.File(path, "r") as h:
        node = h["obs"][column]
        cats = list(node["categories"].asstr()[:])
        codes = node["codes"][:]
    return cats, codes


def _clean_aliases(pairs: list[tuple[str, str]], official: set[str]) -> pd.DataFrame:
    """Keep alias -> symbol pairs whose alias is not an official symbol and points at exactly one symbol."""
    official_lower = {s.lower() for s in official}
    targets: dict[str, set[str]] = defaultdict(set)
    spelled: dict[str, str] = {}
    for alias, symbol in pairs:
        alias = alias.strip()
        if not alias or alias.lower() == symbol.lower() or alias.lower() in official_lower:
            continue
        targets[alias.lower()].add(symbol)
        spelled.setdefault(alias.lower(), alias)
    rows = [(spelled[a], next(iter(s))) for a, s in targets.items() if len(s) == 1]
    return pd.DataFrame(rows, columns=["alias", "symbol"])


def read_mgi_aliases(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df = df[df["Marker Type"].isin(["Gene", "Pseudogene"])]
    pairs = []
    for symbol, syn in zip(df["Marker Symbol"], df["Marker Synonyms (pipe-separated)"]):
        for alias in syn.split("|"):
            pairs.append((alias, symbol))
    return _clean_aliases(pairs, set(df["Marker Symbol"]))


def read_hgnc_aliases(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     usecols=["symbol", "status", "alias_symbol", "prev_symbol"])
    df = df[df["status"] == "Approved"]
    pairs = []
    for symbol, alias_field, prev_field in zip(df["symbol"], df["alias_symbol"], df["prev_symbol"]):
        for alias in (alias_field + "|" + prev_field).split("|"):
            pairs.append((alias, symbol))
    return _clean_aliases(pairs, set(df["symbol"]))


def m39_sources() -> list[Path]:
    """Deposit h5ad files whose var is GENCODE vM39 (identical in both)."""
    out = []
    for key, name in (("immdict", "immdict_mouse_cellbender_counts.h5ad"),
                      ("perturbseq", "geiger_bmdc_e3_genefull_counts.h5ad")):
        d = get_path("deposits", key)
        if d is not None:
            out.append(d / name)
    return out


def main() -> None:
    cfg = load_config()
    ann = cfg["annotations"]
    out = ann_dir()

    mgi = read_mgi_aliases(resolve(ann["mgi_markers"])).assign(species="mouse", source="MGI")
    hgnc = read_hgnc_aliases(resolve(ann["hgnc"])).assign(species="human", source="HGNC")
    aliases = pd.concat([mgi, hgnc])[["species", "alias", "symbol", "source"]]
    aliases.to_csv(out / "aliases_auto.tsv", sep="\t", index=False)
    print(f"aliases: {len(mgi):,} mouse, {len(hgnc):,} human")

    m38 = read_gtf_genes(resolve(ann["mouse_m38_gtf"]))
    v43 = read_gtf_genes(resolve(ann["human_v43_gtf"]))
    tables = [
        m38.loc[~m38.par_y, ["gene_id", "gene_name"]].assign(species="mouse", source="GENCODE vM38"),
        v43.loc[~v43.par_y, ["gene_id", "gene_name"]].assign(species="human", source="GENCODE v43"),
    ]
    # GENCODE vM39 (both mouse deposits) comes from whichever deposit h5ad is available
    for h5ad in m39_sources():
        if h5ad.exists():
            m39 = read_h5ad_var(h5ad)
            tables.append(m39[["gene_id", "gene_name"]].assign(species="mouse", source="GENCODE vM39"))
            break
    id2sym = pd.concat(tables).rename(columns={"gene_name": "symbol"}).drop_duplicates(["gene_id", "symbol"])
    id2sym.to_parquet(out / "id2symbol.parquet", index=False)
    print(f"id2symbol: {len(id2sym):,} rows")


if __name__ == "__main__":
    main()
