"""What each build step needs and what is already built (used by `scbrowser check` and `scbrowser build`)."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from scbrowser.config import app_data, get_path, load_config

DEPOSIT_FILES = {
    "immdict": ["immdict_mouse_cellbender_counts.h5ad", "mouse_edgeR.tar", "human_pseudobulk_all.h5",
                "human_edgeR.tar", "figure_source_tables.tar", "ligand_name_map.tsv", "celltype_lineage_map.tsv"],
    "perturbseq": ["per_target_de_matrices.h5", "per_target_summary.tsv", "genes.tsv",
                   "empirical_fdr_global_calls.tsv", "geiger_bmdc_e3_genefull_counts.h5ad"],
}
DEPOSIT_LABELS = {"immdict": "Immune / Cytokine Dictionary deposit", "perturbseq": "Perturb-seq deposit"}
ZENODO_FILE_URL = "https://zenodo.org/records/{record}/files/{name}?download=1"
STORES = {
    "annotations": "annotations/id2symbol.parquet",
    "cytokine_mouse": "cytokine_mouse/meta.json",
    "cytokine_human": "cytokine_human/meta.json",
    "perturbseq": "perturbseq/ids.parquet",
    "tabula_muris": "tm/embeddings.json",
    "tm_maps": "tm/maps.json",
}
R_INSTALL_HINT = "install R and the Matrix package, e.g. `conda install -c conda-forge r-base r-matrix`"


def deposit_missing(key: str) -> list[Path]:
    """Files of a deposit that are not where the config says (empty list = complete)."""
    d = get_path("deposits", key)
    if d is None or not d.is_dir():
        return [d if d is not None else Path(f"<deposits.{key} not set>")]
    return [d / f for f in DEPOSIT_FILES[key] if not (d / f).exists()]


def zenodo_record(key: str) -> int | None:
    """Zenodo record a deposit is downloaded from (None: not set, so the files must be placed by hand)."""
    return (load_config().get("zenodo") or {}).get(key)


def r_status() -> tuple[bool, str]:
    """(usable, description) for Rscript with the Matrix package, needed to export the Tabula Muris object."""
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript not found; " + R_INSTALL_HINT
    code = 'cat(R.version.string); if (!requireNamespace("Matrix", quietly = TRUE)) quit(status = 3)'
    try:
        res = subprocess.run([rscript, "-e", code], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"Rscript failed: {exc}"
    if res.returncode == 3:
        return False, f"{res.stdout.strip()} without the Matrix package; " + R_INSTALL_HINT
    if res.returncode != 0:
        return False, f"Rscript failed: {res.stderr.strip()[:200]}"
    return True, f"{res.stdout.strip()} with Matrix"


def store_status() -> dict[str, bool]:
    root = app_data()
    return {name: (root / rel).exists() for name, rel in STORES.items()}
