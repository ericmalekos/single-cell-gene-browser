"""Build the Perturb-seq store from the deposit.

    python -m preprocess.build_perturbseq

Writes app_data/perturbseq/ (the app never reads the deposit itself):
    <matrix>.npy           genes x targets: logFC, logCPM, F (float32); PValue, FDR_within_target, FDR_global
                           (float64, which keeps p-values below 1e-45)
    readouts.parquet       gene, calibration_flag, null_rate_P_lt_0_05 (row order of the matrices)
    targets.parquet        target + per_target_summary.tsv columns (column order of the matrices)
    empirical_fdr.tsv      copy of empirical_fdr_global_calls.tsv
    ids.parquet            symbol, current_symbol, gene_ids, readout_idx, target_idx (-1 = role absent)
    global_p_cutoff.json   largest PValue with FDR_global < t, for the volcano cutoff line
"""
from __future__ import annotations

import json
import shutil

import h5py
import numpy as np
import pandas as pd

from preprocess.annotations import read_h5ad_var
from scbrowser.annot import load_aliases
from scbrowser.config import app_data, get_path
from scbrowser.ids import split_ids

FDR_THRESHOLDS = (0.01, 0.05, 0.1)
MATRICES = {"logFC": np.float32, "logCPM": np.float32, "F": np.float32,
            "PValue": np.float64, "FDR_within_target": np.float64, "FDR_global": np.float64}


def main() -> None:
    deposit = get_path("deposits", "perturbseq")
    out = app_data() / "perturbseq"
    out.mkdir(parents=True, exist_ok=True)

    with h5py.File(deposit / "per_target_de_matrices.h5", "r") as h:
        genes = list(h["genes"].asstr()[:])
        targets = list(h["targets"].asstr()[:])
        m = {k: h[k][:].astype(t) for k, t in MATRICES.items()}
        flag = h["gene_calibration_flag"][:].astype(bool)
        null_rate = h["gene_null_rate_P_lt_0_05"][:]

    # README: the Wdr82 fit reproduces the published table (Zeb2os logFC +1.818, FDR within target 4.2e-14)
    g, t = genes.index("Zeb2os"), targets.index("Wdr82")
    lfc, fdr_w = m["logFC"][g, t], m["FDR_within_target"][g, t]
    assert abs(lfc - 1.818) < 5e-4 and abs(fdr_w / 4.2e-14 - 1) < 0.05, (lfc, fdr_w)

    for k, arr in m.items():
        np.save(out / f"{k}.npy", arr)
    pd.DataFrame({"gene": genes, "calibration_flag": flag, "null_rate_P_lt_0_05": null_rate}).to_parquet(
        out / "readouts.parquet", index=False)
    summary = pd.read_csv(deposit / "per_target_summary.tsv", sep="\t").set_index("target").reindex(targets)
    summary.rename_axis("target").reset_index().to_parquet(out / "targets.parquet", index=False)
    shutil.copyfile(deposit / "empirical_fdr_global_calls.tsv", out / "empirical_fdr.tsv")

    readout_ids = pd.read_csv(deposit / "genes.tsv", sep="\t").set_index("gene")["gene_ids"].map(split_ids)
    var = read_h5ad_var(deposit / "geiger_bmdc_e3_genefull_counts.h5ad")
    name_ids = var.groupby("gene_name")["gene_id"].agg(lambda s: sorted(set(s))).to_dict()
    aliases = load_aliases("mouse")

    def current_symbol(s: str) -> str | None:
        return s if s in name_ids else next(iter(aliases.get(s.lower(), [])), None)

    readout_idx = {s: i for i, s in enumerate(genes)}
    target_idx = {s: i for i, s in enumerate(targets)}
    # a target listed under an outdated symbol (March1) whose current symbol is a readout (Marchf1) is one gene
    for t in targets:
        cur = current_symbol(t)
        if t not in readout_idx and cur in readout_idx and cur not in target_idx:
            target_idx[cur] = target_idx.pop(t)
            print(f"  target {t} merged with readout {cur}")
    rows = []
    for s in genes + [x for x in target_idx if x not in readout_idx]:
        current = current_symbol(s)
        ids = readout_ids.get(s) or name_ids.get(current or s, [])
        rows.append((s, current if current != s else None, "|".join(ids),
                     readout_idx.get(s, -1), target_idx.get(s, -1)))
    ids = pd.DataFrame(rows, columns=["symbol", "current_symbol", "gene_ids", "readout_idx", "target_idx"])
    ids.to_parquet(out / "ids.parquet", index=False)

    unresolved = ids[(ids.gene_ids == "")].symbol.tolist()
    renamed = ids[ids.current_symbol.notna()]
    print(f"symbols: {len(ids):,} ({len(genes):,} readouts, {len(targets):,} targets); "
          f"renamed via aliases {len(renamed)}; without Ensembl ID {len(unresolved)}: {', '.join(unresolved[:20])}")

    cutoffs = {}
    for thr in FDR_THRESHOLDS:
        sig = m["FDR_global"] < thr
        cutoffs[str(thr)] = float(m["PValue"][sig].max()) if sig.any() else None
    (out / "global_p_cutoff.json").write_text(json.dumps(cutoffs, indent=2))
    print("FDR_global p cutoffs:", cutoffs)


if __name__ == "__main__":
    main()
