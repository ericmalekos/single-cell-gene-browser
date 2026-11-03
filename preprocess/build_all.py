"""Download the external inputs and build every store under app_data/.

    scbrowser build [--only downloads,annotations,cytokine,perturbseq,tabula_muris,tm_maps]
    (same as python -m preprocess.build_all)

Without --only, steps whose inputs are missing (a deposit, or R for Tabula Muris) are skipped with a note, so
a partial set of datasets still gives a working app. Steps named in --only must be able to run.
Each build runs in its own process, which returns memory between steps (the Tabula Muris export peaks at
~25 GB in R and ~14 GB in Python).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from preprocess.inputs import (DEPOSIT_FILES, DEPOSIT_LABELS, ZENODO_FILE_URL, deposit_missing, r_status,
                               store_status, zenodo_record)
from scbrowser.config import get_path, load_config, resolve

STEPS = ["downloads", "annotations", "cytokine", "perturbseq", "tabula_muris", "tm_maps"]


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def run(*args: str) -> None:
    log("$ " + " ".join(args))
    subprocess.run(list(args), check=True)


def download(url: str, dest: Path, md5: str | None = None) -> None:
    if dest.exists() and not md5:
        log(f"have {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        # curl resumes a partial download (-C -), which matters for the 5.1 GB Zenodo object
        part = dest.with_name(dest.name + ".part")
        run("curl", "-L", "--fail", "--progress-bar", "-C", "-", "--retry", "5", "-o", str(part), url)
        part.rename(dest)
    if md5:
        from preprocess.build_tabula_muris import md5sum
        got = md5sum(dest)
        if got != md5:
            sys.exit(f"md5 mismatch for {dest}: {got} != {md5}; delete it and run again")
        log(f"md5 ok: {dest.name}")


def download_deposit(key: str) -> None:
    """Fetch the deposit files the build needs that are not there yet, md5-checked against the record's MD5SUMS.txt."""
    folder, record = get_path("deposits", key), zenodo_record(key)
    missing = [f for f in DEPOSIT_FILES[key] if not (folder / f).exists()]
    if not missing:
        log(f"have {DEPOSIT_LABELS[key]} ({folder})")
        return
    log(f"{DEPOSIT_LABELS[key]}: downloading {len(missing)} file(s) from Zenodo record {record}")
    sums = folder / "MD5SUMS.txt"
    download(ZENODO_FILE_URL.format(record=record, name="MD5SUMS.txt"), sums)
    md5s = dict(reversed(line.split(maxsplit=1)) for line in sums.read_text().splitlines() if line.strip())
    for name in missing:
        download(ZENODO_FILE_URL.format(record=record, name=name), folder / name, md5=md5s[name])


def plan(requested: set[str] | None) -> tuple[list[str], list[str]]:
    """Steps to run and notes about the ones skipped."""
    steps = [s for s in STEPS if requested is None or s in requested]
    blockers = {}
    for step, key in (("cytokine", "immdict"), ("perturbseq", "perturbseq")):
        missing = deposit_missing(key)
        # the downloads step fetches a missing deposit when its Zenodo record is set
        if missing and not ("downloads" in steps and zenodo_record(key) and get_path("deposits", key)):
            blockers[step] = f"{DEPOSIT_LABELS[key]} incomplete (missing {missing[0]})"
    ok, r_desc = r_status()
    if not ok:
        blockers["tabula_muris"] = r_desc
    if not store_status()["tabula_muris"] and ("tabula_muris" not in steps or "tabula_muris" in blockers):
        blockers["tm_maps"] = "the Tabula Muris store is not built"
    forced = [s for s in steps if requested is not None and s in blockers]
    if forced:
        sys.exit("cannot run " + "; ".join(f"{s}: {blockers[s]}" for s in forced))
    skipped = [f"{s}: {blockers[s]}" for s in steps if s in blockers]
    return [s for s in steps if s not in blockers], skipped


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", help="comma-separated subset of " + ", ".join(STEPS))
    args = ap.parse_args(argv)
    requested = set(args.only.split(",")) if args.only else None
    if requested and requested - set(STEPS):
        sys.exit(f"unknown steps: {', '.join(sorted(requested - set(STEPS)))}")
    steps, skipped = plan(requested)
    for note in skipped:
        log(f"skipping {note}")
    cfg = load_config()
    py = sys.executable

    if "downloads" in steps:
        ann, tm = cfg["annotations"], cfg["tabula_muris"]
        for key in ("mgi_markers", "hgnc", "mouse_m38_gtf", "human_v43_gtf"):
            download(ann[f"{key}_url"], resolve(ann[key]))
        for step, key in (("cytokine", "immdict"), ("perturbseq", "perturbseq")):
            if step in steps:
                download_deposit(key)
        if "tabula_muris" in steps:
            download(tm["global_rds_url"], resolve(tm["global_rds"]), md5=tm["global_rds_md5"])
            if tm.get("immune_rds") and tm.get("immune_rds_url"):
                download(tm["immune_rds_url"], resolve(tm["immune_rds"]), md5=tm.get("immune_rds_md5"))
    if "annotations" in steps:
        run(py, "-m", "preprocess.annotations")
    if "cytokine" in steps:
        run(py, "-m", "preprocess.build_cytokine", "--species", "mouse")
        run(py, "-m", "preprocess.build_cytokine", "--species", "human")
    if "perturbseq" in steps:
        run(py, "-m", "preprocess.build_perturbseq")
    if "tabula_muris" in steps:
        # the downloads step already checked the md5
        run(py, "-m", "preprocess.build_tabula_muris", *(["--skip-md5"] if "downloads" in steps else []))
    if "tm_maps" in steps:
        run(py, "-m", "preprocess.build_tm_maps")
    log("finished: " + ", ".join(steps) + (f"; skipped {len(skipped)} step(s), see above" if skipped else ""))


if __name__ == "__main__":
    main()
