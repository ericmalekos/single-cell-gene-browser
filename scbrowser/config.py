"""Settings: config.yaml defaults merged with config.local.yaml and environment variables.

The config file is $SCBROWSER_CONFIG if set, else config.yaml of the repository (editable install), else
./config.yaml. A config.local.yaml next to it overrides any of its keys. $SCBROWSER_DATA overrides app_data.
Relative paths resolve against the config file's directory.
"""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent


@lru_cache(maxsize=None)
def config_path() -> Path:
    env = os.environ.get("SCBROWSER_CONFIG")
    if env:
        return Path(env).expanduser().resolve()
    for candidate in (ROOT / "config.yaml", Path.cwd() / "config.yaml"):
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError("config.yaml not found: run from the repository or set SCBROWSER_CONFIG")


def _merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


@lru_cache(maxsize=None)
def load_config() -> dict:
    path = config_path()
    cfg = yaml.safe_load(path.read_text()) or {}
    local = path.with_name("config.local.yaml")
    if local.exists():
        cfg = _merge(cfg, yaml.safe_load(local.read_text()) or {})
    return cfg


def local_config_path() -> Path | None:
    local = config_path().with_name("config.local.yaml")
    return local if local.exists() else None


def resolve(p: str | os.PathLike) -> Path:
    """Resolve a config path relative to the config file's directory."""
    p = Path(p).expanduser()
    return p if p.is_absolute() else (config_path().parent / p).resolve()


def get_path(*keys: str) -> Path | None:
    """Path stored under nested keys; None when the setting is empty (e.g. no immune subset)."""
    node = load_config()
    for k in keys:
        node = node[k]
    return None if node in (None, "") else resolve(node)


def app_data() -> Path:
    env = os.environ.get("SCBROWSER_DATA")
    return Path(env).expanduser().resolve() if env else get_path("app_data")


def dataset_info(key: str) -> dict:
    return load_config()["datasets"][key]


def clear_cache() -> None:
    """Forget the loaded settings (tests, or after editing the config files)."""
    config_path.cache_clear()
    load_config.cache_clear()
