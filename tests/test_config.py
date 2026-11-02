"""Config discovery, config.local.yaml overrides and environment variables."""
import pytest

from scbrowser import config


@pytest.fixture
def tmp_config(tmp_path, monkeypatch):
    (tmp_path / "config.yaml").write_text(
        "app_data: app_data\n"
        "deposits:\n  immdict: data/erv_immdict_deposit\n  perturbseq: data/geiger_bmdc_deposit\n"
        "tabula_muris:\n  immune_rds: null\n  check_genes: [Gapdh]\n")
    monkeypatch.setenv("SCBROWSER_CONFIG", str(tmp_path / "config.yaml"))
    monkeypatch.delenv("SCBROWSER_DATA", raising=False)
    config.clear_cache()
    yield tmp_path
    config.clear_cache()


def test_defaults_resolve_relative_to_config(tmp_config):
    assert config.config_path() == tmp_config / "config.yaml"
    assert config.get_path("deposits", "immdict") == tmp_config / "data" / "erv_immdict_deposit"
    assert config.app_data() == tmp_config / "app_data"
    assert config.get_path("tabula_muris", "immune_rds") is None
    assert config.local_config_path() is None


def test_local_overrides_are_merged(tmp_config):
    (tmp_config / "config.local.yaml").write_text(
        "deposits:\n  immdict: /elsewhere/immdict\n"
        "tabula_muris:\n  immune_rds: subset.rds\n")
    config.clear_cache()
    assert config.get_path("deposits", "immdict").as_posix() == "/elsewhere/immdict"
    assert config.get_path("deposits", "perturbseq") == tmp_config / "data" / "geiger_bmdc_deposit"  # kept
    assert config.get_path("tabula_muris", "immune_rds") == tmp_config / "subset.rds"
    assert config.load_config()["tabula_muris"]["check_genes"] == ["Gapdh"]                         # kept


def test_data_dir_from_environment(tmp_config, tmp_path_factory, monkeypatch):
    shared = tmp_path_factory.mktemp("shared")
    monkeypatch.setenv("SCBROWSER_DATA", str(shared))
    assert config.app_data() == shared
