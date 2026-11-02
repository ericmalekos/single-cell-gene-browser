"""Store values against the deposit source tables. Skipped when app_data/ has not been built."""
import tarfile

import numpy as np
import pandas as pd
import pytest

from scbrowser.config import app_data, get_path

pytestmark = pytest.mark.skipif(not (app_data() / "cytokine_mouse").exists(), reason="app_data not built")


def row_of(store, symbol):
    res = store.index.resolve(symbol)
    assert res.found, symbol
    return res.rows[0]


def test_mouse_values_match_panel_table():
    from scbrowser.stores import cytokine_store
    st = cytokine_store("mouse")
    with tarfile.open(get_path("deposits", "immdict") / "mouse_edgeR.tar") as tar:
        panel = pd.read_csv(tar.extractfile("mouse_edgeR/panel/de_B_cell.tsv"), sep="\t")
    d = st.gene_long(row_of(st, "Stat1"))
    for _, ref in panel[panel.gene == "Stat1"].head(20).iterrows():
        got = d[(d.ct == "B cell") & (d.label == ref.cytokine)].iloc[0]
        assert abs(got.logFC - ref.logFC) < 1e-3
        assert abs(got.PValue - ref.PValue) <= 1e-2 * ref.PValue


def test_mouse_pbs_pct_close_to_figure_table():
    from scbrowser.stores import cytokine_store
    st = cytokine_store("mouse")
    with tarfile.open(get_path("deposits", "immdict") / "figure_source_tables.tar") as tar:
        ref = pd.read_csv(tar.extractfile("figure_source_tables/mouse_pct_cells_by_celltype_PBS.tsv"), sep="\t",
                          index_col=0)
    d = st.gene_long(row_of(st, "D17H6S56E-5")).drop_duplicates("ct")
    labels = dict(zip(st.celltypes.display, st.celltypes.label))
    diffs = [abs(r.pct_pbs - ref.loc["D17H6S56E-5", labels[r.ct]]) for r in d.itertuples()]
    assert max(diffs) < 7.5          # README: recomputed values differ by at most 7.3 points


def test_human_overall_expressing_counts_exact():
    from scbrowser.stores import cytokine_store
    st = cytokine_store("human")
    with tarfile.open(get_path("deposits", "immdict") / "figure_source_tables.tar") as tar:
        ref = pd.read_csv(tar.extractfile("figure_source_tables/human_pct_cells.tsv"), sep="\t").set_index("gene")
    for gene in ["TSPAN6", "STAT1", "IRF7", "ERVH48-1"]:
        r = row_of(st, gene)
        total = int(round(float((np.asarray(st.pct[r]) * st.n_cells / 100).sum())))
        assert total == ref.loc[gene, "n_cells_expressing"], gene


def test_untested_genes_have_no_de_values():
    from scbrowser.stores import cytokine_store
    st = cytokine_store("mouse")
    r = int(np.flatnonzero(st.de_row < 0)[0])
    d = st.gene_long(r)
    assert d.logFC.isna().all() and d.pct.notna().any()


@pytest.mark.skipif(not (app_data() / "perturbseq" / "ids.parquet").exists(), reason="Perturb-seq store not built")
def test_perturbseq_wdr82_zeb2os():
    from scbrowser.stores import perturb_store
    st = perturb_store()
    gi, _ = st.roles(row_of(st, "Zeb2os"))
    df = st.readout_table(gi).set_index("label")
    assert abs(df.loc["Wdr82", "logFC"] - 1.818) < 5e-4
    assert abs(df.loc["Wdr82", "FDR_within_target"] / 4.2e-14 - 1) < 0.05
    assert df.PValue.idxmin() == "Wdr82"
    _, ti = st.roles(row_of(st, "Wdr82"))
    assert ti >= 0 and st.within_target_cutoff(ti, 0.05) > 0
    gi, ti = st.roles(row_of(st, "Marchf1"))     # target March1 merged with readout Marchf1
    assert gi >= 0 and st.targets[ti] == "March1"


def gaplinc_stats():
    """The group stats behind Figure_examples/Gaplinc_global.png, if this machine has the Tabula Muris project."""
    immune = get_path("tabula_muris", "immune_rds")
    path = immune.parents[2] / "Figures" / "Gaplinc" / "tissue_violin_group_stats.tsv" if immune else None
    return path if path is not None and path.exists() else None


@pytest.mark.skipif(not (app_data() / "tm" / "cells_global.parquet").exists() or gaplinc_stats() is None,
                    reason="global TM store or the reference group stats not available")
def test_tabula_muris_matches_existing_group_stats():
    from scbrowser.stores import group_stats, tabula_muris_store
    tm = tabula_muris_store()
    cells = tm.embeddings["global"].cells
    v = tm.expression(row_of(tm, "Gm20703"), "global")
    m = (cells.age == "3m").to_numpy()
    mine = group_stats(v[m], cells.tissue[m]).set_index("group")
    ref = pd.read_csv(gaplinc_stats(), sep="\t").set_index("group")
    j = ref[["total_cells", "n_expr_pos"]].join(mine[["cells", "expressing"]])
    # the figures came from an older 3m-only object: the same expressing cells; the current object keeps every
    # cell of it and a few more in some tissues (e.g. 257 more in Spleen), none of which express the gene
    assert (j.expressing == j.n_expr_pos).all()
    assert (j.cells >= j.total_cells).all()


@pytest.mark.skipif(not (app_data() / "tm" / "cells_immune.parquet").exists(), reason="immune TM store not built")
def test_immune_embedding_reuses_global_values():
    from scbrowser.stores import tabula_muris_store
    tm = tabula_muris_store()
    imm = tm.embeddings["immune"]
    assert imm.expr_index is not None and len(imm.cells) == 30411
    r = row_of(tm, "Ptprc")
    glob = tm.expression(r, "global")
    assert np.array_equal(tm.expression(r, "immune"), glob[imm.expr_index])
    assert list(imm.cells.age.cat.categories) == ["3m", "18m", "24m"]
