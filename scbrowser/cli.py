"""Command line: `scbrowser check`, `scbrowser build`, `scbrowser serve`."""
from __future__ import annotations

import argparse
import os
import socket
import subprocess
import sys
from pathlib import Path

from . import __version__
from .config import app_data, config_path, get_path, local_config_path


def _row(status: str, label: str, detail: str = "") -> None:
    print(f"  {status:<9} {label:<38} {detail}")


def cmd_check(args) -> int:
    from preprocess.inputs import DEPOSIT_LABELS, deposit_missing, r_status, store_status

    local = local_config_path()
    print(f"scbrowser {__version__}")
    print(f"config    {config_path()}" + (f" (+ {local.name})" if local else " (no config.local.yaml)"))
    print(f"app_data  {app_data()}\n")
    print("Inputs (needed only to build)")
    for key, label in DEPOSIT_LABELS.items():
        missing = deposit_missing(key)
        if missing:
            _row("MISSING", label, f"{missing[0]}" + (f" (+{len(missing) - 1} more)" if len(missing) > 1 else ""))
        else:
            _row("ok", label, str(get_path("deposits", key)))
    ok, desc = r_status()
    _row("ok" if ok else "MISSING", "R with Matrix (Tabula Muris export)", desc)
    rds = get_path("tabula_muris", "global_rds")
    _row("ok" if rds.exists() else "download", "Tabula Muris object (5.1 GB)", str(rds))
    immune = get_path("tabula_muris", "immune_rds")
    if immune is None:
        _row("optional", "Tabula Muris immune subset", "not configured")
    else:
        _row("ok" if immune.exists() else "MISSING", "Tabula Muris immune subset", str(immune))
    print("\nStores (what the app reads)")
    stores = store_status()
    for name, built in stores.items():
        _row("built" if built else "not built", name)
    if not all(stores.values()):
        print("\nRun `scbrowser build` to create the missing stores (steps without inputs are skipped).")
    return 0


def cmd_build(args) -> int:
    from preprocess.build_all import main as build_all
    build_all(["--only", args.only] if args.only else [])
    return 0


def _local_ip() -> str | None:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("192.0.2.1", 80))   # no packet is sent; this only picks the outgoing interface
            return s.getsockname()[0]
    except OSError:
        return None


def cmd_serve(args, extra: list[str]) -> int:
    from preprocess.inputs import store_status

    built = [k for k, v in store_status().items() if v]
    if not built:
        print("warning: no stores built yet; run `scbrowser build` first", file=sys.stderr)
    app = Path(__file__).with_name("app.py")
    cmd = [sys.executable, "-m", "panel", "serve", str(app), "--port", str(args.port), "--address", args.address]
    origins = list(args.origin)
    if args.address in ("localhost", "127.0.0.1"):
        origins += ["localhost", "127.0.0.1"]
    elif not origins:
        origins += [h for h in ("localhost", socket.getfqdn(), _local_ip()) if h]
    for host in dict.fromkeys(origins):
        cmd += ["--allow-websocket-origin", f"{host}:{args.port}"]
    if args.show:
        cmd.append("--show")
    cmd += extra
    print("starting: " + " ".join(cmd), flush=True)
    if os.name == "posix":
        # become the server process, so its PID (used by the Mac launcher to stop it) and Ctrl+C reach the server
        os.execv(sys.executable, cmd)
    return subprocess.call(cmd)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scbrowser", description="Single-cell gene browser")
    ap.add_argument("--version", action="version", version=f"scbrowser {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("check", help="show the configuration, available inputs and built stores")
    b = sub.add_parser("build", help="download annotations and build the data stores")
    b.add_argument("--only", help="comma-separated: downloads,annotations,cytokine,perturbseq,tabula_muris")
    s = sub.add_parser("serve", help="start the web app; unknown options are passed to `panel serve`")
    s.add_argument("--port", type=int, default=5006)
    s.add_argument("--address", default="localhost", help="0.0.0.0 accepts connections from other machines")
    s.add_argument("--origin", action="append", default=[],
                   help="host name users type in the browser (repeatable; default: this machine's names)")
    s.add_argument("--show", action="store_true", help="open the app in a browser tab")
    args, extra = ap.parse_known_args(argv)
    if extra and args.command != "serve":
        ap.error(f"unrecognized arguments: {' '.join(extra)}")
    if args.command == "check":
        return cmd_check(args)
    if args.command == "build":
        return cmd_build(args)
    return cmd_serve(args, extra)


if __name__ == "__main__":
    sys.exit(main())
