"""Data stores, one per process, opened lazily and shared by every session.

`panel serve` re-runs app.py for each session, so the stores live here (an imported module) rather than in the
app script. Arrays are numpy memmaps, which need no locking and survive forked workers.
"""
from __future__ import annotations

import json
import threading
from dataclasses import dataclass
from functools import lru_cache, wraps
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from .annot import foreign_ids, load_aliases
from .config import app_data
from .genes import GeneIndex
from .ids import split_ids

_lock = threading.RLock()


def _singleton(fn):
    cached = lru_cache(maxsize=None)(fn)

    @wraps(fn)
    def wrapper(*args):
        with _lock:
            return cached(*args)

    wrapper.cache_clear = cached.cache_clear
    return wrapper


# -- Tabula Muris ------------------------------------------------------------------------------------------------
@dataclass
class Embedding:
    key: str
    label: str
    cells: pd.DataFrame
    expr: str                      # expression store directory
    expr_index: np.ndarray | None  # cell -> column of the expression store (None = identity)
    map: Path | None = None        # labelled reference map (PNG), if built
    map_colored_by: str = ""
    map_size: tuple[int, int] = (0, 0)   # pixels


class TabulaMurisStore:
    species = "mouse"

    def __init__(self):
        root = app_data() / "tm"
        self.root = root
        self.genes = pd.read_parquet(root / "genes.parquet")
        maps = json.loads((root / "maps.json").read_text()) if (root / "maps.json").exists() else {}
        self.embeddings: dict[str, Embedding] = {}
        for e in json.loads((root / "embeddings.json").read_text()):
            cells = pd.read_parquet(root / e["cells"])
            idx = cells.pop("expr_index").to_numpy(np.int64) if "expr_index" in cells else None
            emb = Embedding(e["key"], e["label"], cells, e["expr"], idx)
            m = maps.get(e["key"])
            if m and (root / m["file"]).exists():
                with Image.open(root / m["file"]) as im:    # reads the header only
                    emb.map, emb.map_colored_by, emb.map_size = root / m["file"], m["colored_by"], im.size
            self.embeddings[e["key"]] = emb
        self._expr: dict[str, tuple] = {}
        self.index = GeneIndex(self.genes.symbol, [split_ids(x) for x in self.genes.gene_ids], species="mouse",
                               hidden=self.genes.is_spike.tolist(), aliases=load_aliases("mouse"),
                               foreign_ids=foreign_ids)

    def _arrays(self, name: str):
        if name not in self._expr:
            d = self.root / name
            shape = json.loads((d / "shape.json").read_text())
            self._expr[name] = (np.load(d / "indptr.npy", mmap_mode="r"), np.load(d / "indices.npy", mmap_mode="r"),
                                np.load(d / "data.npy", mmap_mode="r"), shape["n_cells"])
        return self._expr[name]

    def expression(self, row: int, embedding: str) -> np.ndarray:
        emb = self.embeddings[embedding]
        indptr, indices, data, n_cells = self._arrays(emb.expr)
        s, e = int(indptr[row]), int(indptr[row + 1])
        v = np.zeros(n_cells, np.float32)
        v[indices[s:e]] = data[s:e]
        return v if emb.expr_index is None else v[emb.expr_index]


def group_stats(values: np.ndarray, groups: pd.Series) -> pd.DataFrame:
    """Per-group summary as in the *_group_stats.tsv files, plus % expressing."""
    df = pd.DataFrame({"group": groups.astype(str).to_numpy(), "v": values})
    df = df[df.group != ""]
    g = df.groupby("group").v
    out = pd.DataFrame({
        "cells": g.size(),
        "expressing": g.apply(lambda s: int((s > 0).sum())),
        "mean": g.mean(),
        "median": g.median(),
        "q3": g.quantile(0.75),
        "max": g.max(),
    })
    out.insert(2, "pct_expressing", 100 * out.expressing / out.cells)
    return out.reset_index().sort_values(["mean", "pct_expressing"], ascending=False, ignore_index=True)


# -- cytokine dictionaries ---------------------------------------------------------------------------------------
class CytokineStore:
    def __init__(self, species: str):
        root = app_data() / f"cytokine_{species}"
        self.species = species
        self.genes = pd.read_parquet(root / "genes.parquet")
        self.celltypes = pd.read_parquet(root / "celltypes.parquet")
        self.cytokines = pd.read_parquet(root / "cytokines.parquet")
        self.meta = json.loads((root / "meta.json").read_text())
        self.logfc = np.load(root / "logFC.npy", mmap_mode="r")
        self.pval = np.load(root / "PValue.npy", mmap_mode="r")
        self.fdr = np.load(root / "FDR.npy", mmap_mode="r")
        self.pct = np.load(root / "pct_expr.npy", mmap_mode="r")
        self.n_cells = np.load(root / "n_cells.npy")
        self.de_row = self.genes.de_row.to_numpy()
        variant_of = self.genes["variant_of"].tolist() if "variant_of" in self.genes else None
        self.index = GeneIndex(self.genes.symbol, [split_ids(x) for x in self.genes.gene_ids], species=species,
                               variant_of=variant_of, aliases=load_aliases(species), foreign_ids=foreign_ids)

    def gene_long(self, row: int) -> pd.DataFrame:
        """One row per (cell type, cytokine) for one gene."""
        n_ct, n_cyt = len(self.celltypes), len(self.cytokines)
        r = int(self.de_row[row])
        pct = np.asarray(self.pct[row], np.float64)               # (n_ct, n_cyt + 1)
        if r >= 0:
            lfc, p, q = (np.asarray(a[r], np.float64) for a in (self.logfc, self.pval, self.fdr))
        else:
            lfc = p = q = np.full((n_ct, n_cyt), np.nan)
        ct = np.repeat(self.celltypes.display.to_numpy(), n_cyt)
        cyt = np.tile(self.cytokines.display.to_numpy(), n_ct)
        df = pd.DataFrame({
            "ct": ct, "cyt": cyt, "logFC": lfc.ravel(), "PValue": p.ravel(), "FDR": q.ravel(),
            "pct": pct[:, :n_cyt].ravel(), "pct_pbs": np.repeat(pct[:, n_cyt], n_cyt),
            "n_cells": self.n_cells[:, :n_cyt].ravel(), "n_cells_pbs": np.repeat(self.n_cells[:, n_cyt], n_cyt),
        })
        fam = self.cytokines.set_index("display")
        for col in ("DEG_plot_family", "label"):
            df[col] = df.cyt.map(fam[col])
        df["lineage"] = df.ct.map(self.celltypes.set_index("display").lineage)
        return df


# -- Perturb-seq -------------------------------------------------------------------------------------------------
class PerturbStore:
    species = "mouse"
    MATRICES = ("logFC", "logCPM", "F", "PValue", "FDR_within_target", "FDR_global")

    def __init__(self):
        root = app_data() / "perturbseq"
        self.m = {k: np.load(root / f"{k}.npy", mmap_mode="r") for k in self.MATRICES}   # genes x targets
        readouts = pd.read_parquet(root / "readouts.parquet")
        self.readouts = readouts.gene.tolist()
        self.flag = readouts.calibration_flag.to_numpy(bool)
        self.null_rate = readouts.null_rate_P_lt_0_05.to_numpy()
        self.summary = pd.read_parquet(root / "targets.parquet").set_index("target")
        self.targets = self.summary.index.tolist()
        self.empirical = pd.read_csv(root / "empirical_fdr.tsv", sep="\t")
        self.cutoffs = {float(k): v for k, v in json.loads((root / "global_p_cutoff.json").read_text()).items()}
        self.ids = pd.read_parquet(root / "ids.parquet")
        self.index = GeneIndex(self.ids.symbol, [split_ids(x) for x in self.ids.gene_ids], species="mouse",
                               current=self.ids.current_symbol.tolist(), aliases=load_aliases("mouse"),
                               foreign_ids=foreign_ids)
        tgt = self.ids[self.ids.target_idx >= 0]
        self.target_rows = np.full(len(self.targets), -1, np.int64)   # KO target -> its row in ids
        self.target_rows[tgt.target_idx.to_numpy()] = tgt.index.to_numpy()

    def targets_overview(self) -> pd.DataFrame:
        """One row per KO target, A-Z: current symbol when renamed, cells, and how many genes it changes."""
        rows = self.ids.iloc[self.target_rows]
        listed = np.asarray(self.targets, dtype=object)
        df = self.summary.reset_index()
        df.insert(1, "current_symbol", np.where(rows.symbol.to_numpy() != listed, rows.symbol.to_numpy(),
                                                rows.current_symbol.fillna("").to_numpy()))
        df["readout"] = rows.readout_idx.to_numpy() >= 0
        df["ti"] = np.arange(len(df))
        counts = [c for c in df.columns if c.startswith(("n_up", "n_down", "n_channels", "n_cells"))]
        df[counts] = df[counts].astype(np.int64)
        return df.sort_values("target", key=lambda t: t.str.lower(), ignore_index=True)

    def roles(self, row: int) -> tuple[int, int]:
        r = self.ids.iloc[row]
        return int(r.readout_idx), int(r.target_idx)

    def readout_table(self, gi: int) -> pd.DataFrame:
        """Every KO's effect on readout gene gi."""
        df = pd.DataFrame({"label": self.targets, **{k: np.asarray(v[gi, :]) for k, v in self.m.items()}})
        for col in ("n_channels", "n_cells_KO", "n_cells_NTC"):
            df[col] = self.summary[col].to_numpy()
        return df

    def target_table(self, ti: int) -> pd.DataFrame:
        """Effect of knocking out target ti on every readout gene."""
        df = pd.DataFrame({"label": self.readouts, **{k: np.asarray(v[:, ti]) for k, v in self.m.items()}})
        df["calibration_flag"] = self.flag
        df["null_rate_P_lt_0_05"] = self.null_rate
        return df

    def within_target_cutoff(self, ti: int, thr: float) -> float | None:
        sig = np.asarray(self.m["FDR_within_target"][:, ti]) < thr
        return float(np.asarray(self.m["PValue"][:, ti])[sig].max()) if sig.any() else None


@_singleton
def tabula_muris_store() -> TabulaMurisStore:
    return TabulaMurisStore()


@_singleton
def cytokine_store(species: str) -> CytokineStore:
    return CytokineStore(species)


@_singleton
def perturb_store() -> PerturbStore:
    return PerturbStore()
