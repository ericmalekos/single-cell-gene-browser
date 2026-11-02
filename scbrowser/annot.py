"""Readers for the generated annotation tables (aliases, cross-dataset Ensembl ID -> symbol)."""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

import pandas as pd

from .config import app_data, get_path


@lru_cache(maxsize=None)
def load_aliases(species: str) -> dict[str, list[str]]:
    """lower-case alias -> [symbol]; manual aliases.tsv entries override the generated MGI/HGNC table."""
    out: dict[str, list[str]] = {}
    sources = [app_data() / "annotations" / "aliases_auto.tsv", get_path("aliases")]
    for path in sources:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        df = df[df.species == species]
        for alias, symbol in zip(df.alias, df.symbol):
            out[alias.lower()] = [symbol]
    return out


@lru_cache(maxsize=None)
def _id_table() -> dict[str, list[tuple[str, str]]]:
    path = app_data() / "annotations" / "id2symbol.parquet"
    table: dict[str, list[tuple[str, str]]] = defaultdict(list)
    if path.exists():
        df = pd.read_parquet(path, columns=["gene_id", "symbol", "species"])
        for gid, sym, sp in zip(df.gene_id, df.symbol, df.species):
            if (sym, sp) not in table[gid]:
                table[gid].append((sym, sp))
    return dict(table)


def foreign_ids(gene_id: str) -> list[tuple[str, str]]:
    """[(symbol, species)] for an unversioned Ensembl gene ID from GENCODE vM38, vM39 or v43."""
    return _id_table().get(gene_id, [])
