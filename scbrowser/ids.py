"""Ensembl gene ID helpers shared by preprocessing and the app."""
from __future__ import annotations

import re

# ENSG00000141510, ENSMUSG00000026104.16, ENSG00000002586.20_PAR_Y
_ENS_RE = re.compile(r"^(ENS[A-Z]{0,6}G\d{11})(?:\.\d+)?(?:_PAR_Y)?$", re.IGNORECASE)


def is_ensembl(text: str) -> bool:
    return bool(_ENS_RE.match(text.strip()))


def strip_ens(text: str) -> str:
    """Drop the version and any _PAR_Y suffix; return the input unchanged if it is not an Ensembl gene ID."""
    m = _ENS_RE.match(text.strip())
    return m.group(1).upper() if m else text.strip()


def split_ids(value) -> list[str]:
    """Split a '|' or ',' separated ID field into unversioned IDs."""
    if value is None or (isinstance(value, float) and value != value):
        return []
    parts = re.split(r"[|,]", str(value))
    return [strip_ens(p) for p in parts if p.strip()]
