"""Cytokine dictionary view (mouse or human): dot plot of responses to every cytokine in every cell type."""
from __future__ import annotations

import numpy as np
import pandas as pd
import panel as pn

from ..config import dataset_info
from ..plots.dotplot import dotplot_figure
from ..stores import cytokine_store
from .base import FMT_FC, FMT_NUM, FMT_P, FMT_PCT, DatasetView, card_html, esc, fmt_p

FAMILY_COLUMNS = {"Manuscript families": "DEG_plot_family", "Single family": "single_family",
                  "Family level 1": "family_level_1", "Family level 2": "family_level_2"}
SUBSETS = ["All", "Significant", "Families", "Custom"]


class CytokineView(DatasetView):
    def __init__(self, species: str):
        self.key = f"{species}_cytokine"
        super().__init__(dataset_info(self.key))
        self.species = species
        self.store = cytokine_store(species)
        cyts = self.store.cytokines
        self.family_col = pn.widgets.Select(label="Group cytokines by", options=FAMILY_COLUMNS,
                                            value="DEG_plot_family", width=190)
        self.subset = pn.widgets.RadioButtonGroup(label="Cytokines", options=SUBSETS, value="All",
                                                  color="primary", variant="outline")
        self.families = pn.widgets.MultiChoice(label="Families", options=self._family_options(), visible=False,
                                               width=420)
        self.custom = pn.widgets.MultiChoice(label="Cytokines", options=list(cyts.display), visible=False, width=420)
        self.fdr = pn.widgets.Select(label="FDR threshold", options=[0.01, 0.05, 0.1], value=0.05, width=120)
        self.clim = pn.widgets.Select(label="Color limit |log2FC|", value="Auto", width=150,
                                      options=["Auto", "1", "2", "3", "4", "6"])
        self.sig_only = pn.widgets.Checkbox(label="Significant rows only", value=False)
        self.plot = pn.pane.Bokeh(sizing_mode="fixed")
        self.message = pn.pane.HTML("")
        self._long: pd.DataFrame | None = None
        for w in (self.subset, self.families, self.custom, self.fdr, self.clim):
            w.param.watch(lambda e: self._redraw(), "value")
        self.family_col.param.watch(self._family_col_changed, "value")
        self.sig_only.param.watch(lambda e: self._update_table(), "value")

    def _family_options(self) -> list[str]:
        col = self.family_col.value if hasattr(self, "family_col") else "DEG_plot_family"
        cyts = self.store.cytokines
        return list(cyts.groupby(col)["order"].min().sort_values().index)

    def _family_col_changed(self, event) -> None:
        self.families.value = []
        self.families.options = self._family_options()
        self._redraw()

    def panel(self):
        controls = pn.Row(pn.Column(pn.widgets.StaticText(value="Cytokines"), self.subset), self.family_col,
                          self.fdr, self.clim, sizing_mode="stretch_width")
        plot_box = pn.Column(self.plot, styles={"overflow-x": "auto"}, sizing_mode="stretch_width")
        return pn.Column(self.header, self.card, controls, pn.Row(self.families, self.custom), self.message,
                         plot_box, pn.Row(self.sig_only, pn.Spacer(), self.download, sizing_mode="stretch_width"),
                         self.table, *self.methods, sizing_mode="stretch_width")

    def show(self, row, resolution=None):
        self.row, self.resolution = row, resolution
        self._long = self.store.gene_long(row)
        self._redraw()

    def clear_plots(self) -> None:
        self._long = None
        self.plot.object = None
        self.message.object = ""

    def _selected_cytokines(self) -> pd.DataFrame:
        cyts, mode = self.store.cytokines, self.subset.value
        self.families.visible = mode == "Families"
        self.custom.visible = mode == "Custom"
        if mode == "Significant":
            sig = self._long.loc[self._long.FDR < self.fdr.value, "cyt"].unique()
            return cyts[cyts.display.isin(sig)]
        if mode == "Families":
            return cyts[cyts[self.family_col.value].isin(self.families.value)]
        if mode == "Custom":
            return cyts[cyts.display.isin(self.custom.value)]
        return cyts

    def _redraw(self) -> None:
        if self._long is None:
            return
        sym = self.symbol(self.row)
        chosen = self._selected_cytokines()
        if chosen.empty:
            hint = {"Significant": f"No cytokine reaches FDR < {self.fdr.value:g} in any cell type.",
                    "Families": "Choose one or more families.", "Custom": "Choose one or more cytokines."}
            self.message.object = f"<div class='note'>{hint.get(self.subset.value, '')}</div>"
            self.plot.object = None
        else:
            self.message.object = ""
            clim = None if self.clim.value == "Auto" else float(self.clim.value)
            title = f"{sym}: response to cytokines vs PBS ({self.info['title']})"
            self.plot.object = dotplot_figure(self._long, list(self.store.celltypes.display), chosen,
                                              family_col=self.family_col.value, fdr_thr=self.fdr.value, clim=clim,
                                              title=title)
        self._update_card(sym)
        self._update_table()

    def _update_card(self, sym: str) -> None:
        d, thr = self._long, self.fdr.value
        tested_ct = d.groupby("ct").logFC.apply(lambda s: s.notna().any())
        sig = d[d.FDR < thr]
        stats = [("cell types tested", f"{int(tested_ct.sum())} / {len(tested_ct)}"),
                 (f"up, FDR < {thr:g}", f"{int((sig.logFC > 0).sum())}"),
                 (f"down, FDR < {thr:g}", f"{int((sig.logFC < 0).sum())}"),
                 ("max % expressing (PBS)", f"{np.nanmax(d.pct_pbs):.1f}%" if d.pct_pbs.notna().any() else "NA")]
        notes, warnings = [], []
        if len(sig):
            top = sig.reindex(sig.logFC.abs().sort_values(ascending=False).index).head(4)
            notes.append("Strongest significant responses: " + "; ".join(
                f"{esc(r.cyt)} in {esc(r.ct)} (log2FC {r.logFC:+.2f}, FDR {fmt_p(r.FDR)})" for r in top.itertuples()))
        if not tested_ct.any():
            warnings.append("Not tested in any cell type (did not pass filterByExpr). Dot sizes still show the "
                            "% of cells expressing the gene.")
        notes.append("FDR: Benjamini-Hochberg over all genes within each cell type x cytokine (the deposit's "
                     "genome-wide tables). Manuscript figures adjusted over smaller test sets, so their outlines can "
                     "differ.")
        if self.resolution is not None and self.resolution.note:
            notes.append(esc(self.resolution.note))
        self.card.object = card_html(sym, self.ids_text(self.row), stats, notes, warnings)

    def _update_table(self) -> None:
        if self._long is None:
            return
        d = self._long
        if self.sig_only.value:
            d = d[d.FDR < self.fdr.value]
        table = pd.DataFrame({
            "cell_type": d.ct, "lineage": d.lineage, "cytokine": d.cyt, "family": d.DEG_plot_family,
            "log2FC": d.logFC, "PValue": d.PValue, "FDR": d.FDR, "pct_expressing": d.pct,
            "pct_expressing_PBS": d.pct_pbs, "cells": d.n_cells, "cells_PBS": d.n_cells_pbs,
        }).sort_values("PValue", na_position="last", ignore_index=True)
        self.set_table(table, f"{self.symbol(self.row)}_{self.key}.csv",
                       {"log2FC": FMT_FC, "PValue": FMT_P, "FDR": FMT_P, "pct_expressing": FMT_PCT,
                        "pct_expressing_PBS": FMT_PCT, "cells": FMT_NUM, "cells_PBS": FMT_NUM})
