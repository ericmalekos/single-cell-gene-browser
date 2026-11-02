"""Resolve a user query (symbol, alias or Ensembl ID) to rows of one dataset."""
from __future__ import annotations

import bisect
import difflib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Mapping, Sequence

from .ids import is_ensembl, strip_ens

CROSS_SPECIES_NOTE = "matched by symbol across species; this is not an orthology mapping"


@dataclass
class Resolution:
    query: str
    rows: list[int] = field(default_factory=list)   # best match first; empty when not found
    how: str = ""
    note: str = ""
    suggestions: list[str] = field(default_factory=list)

    @property
    def found(self) -> bool:
        return bool(self.rows)


def _clean(values, n: int) -> list[str | None]:
    """Optional per-row strings; NaN/empty (as read back from parquet) become None."""
    if values is None:
        return [None] * n
    return [v if isinstance(v, str) and v else None for v in values]


class GeneIndex:
    """Lookup tables for one dataset's genes.

    symbols     display symbol per row (duplicates allowed)
    ids         unversioned Ensembl IDs per row
    current     official symbol when the row's symbol is outdated (Perturb-seq 'March1' -> 'Marchf1')
    variant_of  base symbol when the row is an anndata make_unique variant ('MATR3-1' -> 'MATR3')
    hidden      rows kept searchable but left out of completions (ERCC spike-ins)
    aliases     lower-case alias -> official symbols
    foreign_ids Ensembl ID -> [(symbol, species)] from other annotations, for IDs the dataset lacks
    """

    def __init__(
        self,
        symbols: Sequence[str],
        ids: Sequence[Sequence[str]] | None = None,
        *,
        species: str = "",
        current: Sequence[str | None] | None = None,
        variant_of: Sequence[str | None] | None = None,
        hidden: Sequence[bool] | None = None,
        aliases: Mapping[str, Sequence[str]] | None = None,
        foreign_ids: Callable[[str], Sequence[tuple[str, str]]] | None = None,
    ):
        self.symbols = [str(s) for s in symbols]
        n = len(self.symbols)
        self.ids = [list(x) for x in ids] if ids is not None else [[] for _ in range(n)]
        self.species = species
        self.current = _clean(current, n)
        self.variant_of = _clean(variant_of, n)
        self.hidden = list(hidden) if hidden is not None else [False] * n
        self._aliases = {k.lower(): list(v) for k, v in (aliases or {}).items()}
        self._foreign = foreign_ids

        self._exact: dict[str, list[int]] = defaultdict(list)
        self._lower: dict[str, list[int]] = defaultdict(list)
        self._by_id: dict[str, list[int]] = defaultdict(list)
        self._variants: dict[str, list[int]] = defaultdict(list)   # lower base symbol -> variant rows
        self._by_current: dict[str, list[int]] = defaultdict(list)  # lower official symbol -> outdated rows
        for i, s in enumerate(self.symbols):
            self._exact[s].append(i)
            self._lower[s.lower()].append(i)
            for g in self.ids[i]:
                self._by_id[strip_ens(g)].append(i)
            if self.variant_of[i]:
                self._variants[self.variant_of[i].lower()].append(i)
            cur = self.current[i]
            if cur and cur.lower() != s.lower():
                self._by_current[cur.lower()].append(i)

        visible = sorted({s for s, h in zip(self.symbols, self.hidden) if not h}, key=str.lower)
        self._completions = visible
        self._sorted_lower = [s.lower() for s in visible]

    def __len__(self) -> int:
        return len(self.symbols)

    def completions(self) -> list[str]:
        return self._completions

    def label(self, row: int) -> str:
        ids = ", ".join(self.ids[row])
        return f"{self.symbols[row]} ({ids})" if ids else self.symbols[row]

    # -- lookup --------------------------------------------------------------------------------------------------
    def _with_variants(self, rows: list[int]) -> list[int]:
        out = list(rows)
        for r in rows:
            base = (self.variant_of[r] or self.symbols[r]).lower()
            for extra in self._lower.get(base, []) + self._variants.get(base, []):
                if extra not in out:
                    out.append(extra)
        return out

    def _symbol_rows(self, symbol: str) -> tuple[list[int], str]:
        if symbol in self._exact:
            return self._with_variants(self._exact[symbol]), "symbol"
        low = symbol.lower()
        if low in self._lower:
            return self._with_variants(self._lower[low]), "symbol (case-insensitive)"
        if low in self._by_current:
            rows = self._by_current[low]
            listed = ", ".join(sorted({self.symbols[r] for r in rows}))
            return list(rows), f"current symbol; listed in this dataset as {listed}"
        return [], ""

    def resolve(self, query: str) -> Resolution:
        q = (query or "").strip()
        if not q:
            return Resolution(q)

        if is_ensembl(q):
            gid = strip_ens(q)
            if gid in self._by_id:
                return Resolution(q, self._with_variants(self._by_id[gid]), "Ensembl ID")
            for symbol, species in (self._foreign(gid) if self._foreign else []):
                rows, how = self._symbol_rows(symbol)
                if rows:
                    note = CROSS_SPECIES_NOTE if species and self.species and species != self.species else ""
                    return Resolution(q, rows, f"{gid} is {symbol} ({species}); {how}", note)
            return Resolution(q, how="Ensembl ID not in this dataset")

        rows, how = self._symbol_rows(q)
        if rows:
            note = CROSS_SPECIES_NOTE if self._other_species_casing(q, self.symbols[rows[0]]) else ""
            return Resolution(q, rows, how, note)

        for target in self._aliases.get(q.lower(), []):
            rows, how = self._symbol_rows(target)
            if rows:
                return Resolution(q, rows, f"alias of {target}")

        return Resolution(q, suggestions=self.suggest(q))

    def _other_species_casing(self, query: str, matched: str) -> bool:
        """A human-style query (STAT1) hitting a mouse symbol (Stat1), or the reverse."""
        if query == matched or not any(c.isalpha() for c in query):
            return False
        if self.species == "mouse":
            return query.isupper() and not matched.isupper()
        if self.species == "human":
            return not query.isupper() and matched.isupper()
        return False

    def suggest(self, query: str, n: int = 8) -> list[str]:
        low = query.strip().lower()
        if not low:
            return []
        out: list[str] = []
        i = bisect.bisect_left(self._sorted_lower, low)
        while i < len(self._sorted_lower) and self._sorted_lower[i].startswith(low) and len(out) < n:
            out.append(self._completions[i])
            i += 1
        if len(out) < n:
            # difflib over the whole list is slow; compare against symbols sharing the first character
            lo = bisect.bisect_left(self._sorted_lower, low[0])
            hi = bisect.bisect_left(self._sorted_lower, chr(ord(low[0]) + 1))
            pool = {self._sorted_lower[j]: self._completions[j] for j in range(lo, hi)}
            for m in difflib.get_close_matches(low, list(pool), n=n, cutoff=0.75):
                if pool[m] not in out and len(out) < n:
                    out.append(pool[m])
        return out
