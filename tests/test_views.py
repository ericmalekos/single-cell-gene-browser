"""Headless smoke tests: build every view, show known genes, and build the Bokeh models."""
import numpy as np
import pytest
from bokeh.document import Document

from scbrowser.config import app_data

pn = pytest.importorskip("panel")
pytestmark = pytest.mark.skipif(not (app_data() / "cytokine_mouse").exists(), reason="app_data not built")


def show(view, symbol):
    res = view.index.resolve(symbol)
    assert res.found, symbol
    view.show(res.rows[0], res)
    return res


def render(view):
    view.panel().get_root(Document())   # raises on invalid Bokeh property values


@pytest.mark.parametrize("species,genes", [("mouse", ["Synir", "Stat1", "Xist"]), ("human", ["STAT1", "MATR3"])])
def test_cytokine_view(species, genes):
    from scbrowser.views.cytokine import CytokineView
    v = CytokineView(species)
    for g in genes:
        show(v, g)
        assert len(v.table.value) == len(v.store.celltypes) * len(v.store.cytokines)
        assert v.plot.object is not None
    for mode in ("Significant", "Custom", "All"):
        v.subset.value = mode
    v.family_col.value = "family_level_2"
    render(v)


@pytest.mark.skipif(not (app_data() / "perturbseq" / "ids.parquet").exists(), reason="Perturb-seq store not built")
def test_perturbseq_view_modes():
    from scbrowser.views.perturbseq import TARGET, PerturbView
    v = PerturbView()
    show(v, "Zeb2os")
    assert v.mode.options == ["Gene as readout"] and len(v.table.value) == 1130
    show(v, "Rigi")                                  # KO-only target listed as Ddx58
    assert v.mode.value == TARGET and v.fdr_col.value == "FDR_within_target" and len(v.table.value) == 2300
    show(v, "Wdr82")
    v.mode.value = TARGET
    assert v.fdr_col.value == "FDR_within_target"
    render(v)


@pytest.mark.skipif(not (app_data() / "perturbseq" / "ids.parquet").exists(), reason="Perturb-seq store not built")
def test_perturbseq_view_button_follows_new_options():
    from scbrowser.views.perturbseq import TARGET, PerturbView
    v = PerturbView()
    root = v.panel().get_root(Document())
    button = v.mode._models[root.ref["id"]][0]
    show(v, "Rigi")                                  # KO only: options [KO]
    show(v, "Wdr82")                                 # both roles, KO view kept: options [readout, KO]
    assert v.mode.value == TARGET and button.labels[button.active] == TARGET


@pytest.mark.skipif(not (app_data() / "perturbseq" / "ids.parquet").exists(), reason="Perturb-seq store not built")
def test_perturbseq_ko_target_list():
    from scbrowser.views.perturbseq import READOUT, TARGET, PerturbView
    v = PerturbView()
    assert len(v.target_search.options) == len(v.targets_table.value) == 1130
    show(v, "Zeb2os")                                # not knocked out: nothing selected in the KO list
    assert v.target_search.value == "" and v.targets_table.selection == []
    show(v, "Wdr82")                                 # readout and target: the readout view comes first
    assert v.mode.value == READOUT and v.target_search.value == ""
    v.target_search.value = "Wdr82"                  # picking the gene already shown switches to its KO view
    assert v.mode.value == TARGET and v.overview.target.iloc[v.targets_table.selection[0]] == "Wdr82"
    v.mode.value = READOUT                           # the box shows the KO being viewed, so it empties again
    assert v.target_search.value == "" and v.targets_table.selection == []

    picked = []
    v.select_gene = lambda symbol: (picked.append(symbol), show(v, symbol))
    v.target_search.value = "Ddx58 (Rigi)"           # renamed target, listed with its current symbol
    assert picked == ["Ddx58"] and v.store.targets[v._roles[1]] == "Ddx58" and v.mode.value == TARGET
    show(v, "Zeb2os")
    pos = int(np.flatnonzero(v.overview.target == "Rfwd2")[0])   # merged into the readout Cop1
    v._pick_target(int(v.overview.ti.iloc[pos]))
    assert picked[-1] == "Cop1" and v.mode.value == TARGET and v.target_search.value == "Rfwd2 (Cop1)"
    render(v)


@pytest.mark.skipif(not (app_data() / "tm").exists(), reason="TM store not built")
def test_tabula_muris_view():
    from scbrowser.views.tabula_muris import TabulaMurisView
    v = TabulaMurisView()
    show(v, "Gm20703")
    emb = v.store.embeddings[v.embedding.value]
    assert v.map_box.visible == (emb.map is not None)
    if emb.map is not None:
        assert v.map.object == str(emb.map) and emb.map_colored_by.lower() in v.map_note.object
    n_fg = len(v.umap.fg.data["x"])
    assert 0 < n_fg < len(v.store.embeddings[v.embedding.value].cells)
    v.age.value = "3m"
    assert len(v.umap.bg.data["x"]) < len(v.store.embeddings[v.embedding.value].cells)
    if "immune" in v.store.embeddings:
        v.age.value = "All"
        v.embedding.value = "immune"
        assert len(v.umap.bg.data["x"]) == len(v.store.embeddings["immune"].cells)
        imm = v.store.embeddings["immune"]
        assert v.map.object == (str(imm.map) if imm.map else None)
    v.group_by.value = "tissue"
    assert "Tissue" in v.table.value.columns
    render(v)


def test_every_dataset_has_a_methods_section():
    from scbrowser.ui import DATASETS
    from scbrowser.views.base import METHODS_DIR, methods_card
    for key in DATASETS:
        assert (METHODS_DIR / f"{key}.md").exists(), key
        (card,) = methods_card(key)
        md = card.objects[0].object
        assert card.collapsed and "#### References" in md and 'target="_blank"' in md and "](http" not in md
    assert methods_card("no-such-dataset") == []


def test_app_shell_builds():
    from scbrowser.ui import BrowserApp
    app = BrowserApp()
    app._activate()
    view = app.views[app.active_key]
    assert view.row is not None                     # default gene resolved in the first tab
    app.search.value = "definitely-not-a-gene"
    assert not app.chooser.visible and view.row is None
    assert len(view.umap.fg.data["x"]) == 0 and view.violin_cell.object is None   # old gene's figures removed


@pytest.mark.skipif(not (app_data() / "perturbseq" / "ids.parquet").exists(), reason="Perturb-seq store not built")
def test_ko_target_pick_updates_the_whole_app():
    from scbrowser.ui import DATASETS, BrowserApp
    from scbrowser.views.perturbseq import TARGET
    app = BrowserApp()
    app.tabs.active = DATASETS.index("perturbseq")
    view = app.views["perturbseq"]
    view.target_search.value = "Wdr82"
    assert app.search.value == "Wdr82" and view.mode.value == TARGET
    assert "Wdr82" in app.status.object and app.stale == set(DATASETS) - {"perturbseq"}
