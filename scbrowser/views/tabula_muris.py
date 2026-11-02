"""Tabula Muris view: UMAP colored by expression beside a labelled map, top-20 violins, per-group statistics."""
from __future__ import annotations

import numpy as np
import pandas as pd
import panel as pn

from ..config import dataset_info
from ..plots.umap import UmapPlot
from ..plots.violin import violin_figure
from ..stores import group_stats, tabula_muris_store
from .base import FMT_FC, FMT_NUM, FMT_PCT, DatasetView, card_html, esc

GROUPS = {"Cell ontology class": "cell_type", "Tissue": "tissue", "Free annotation": "free_annotation",
          "Cluster": "cluster", "Age": "age"}
GROUP_LABEL = {v: k for k, v in GROUPS.items()}
LABELS = {"None": "", "Tissue": "tissue", "Cell ontology class": "cell_type"}
MAP_WIDTH, MAP_HEIGHT = 720, 500   # fits beside the 760-wide UMAP on a 1920-wide screen


class TabulaMurisView(DatasetView):
    key = "tabula_muris"

    def __init__(self):
        super().__init__(dataset_info(self.key))
        self.store = tabula_muris_store()
        embs = self.store.embeddings
        self.embedding = pn.widgets.Select(label="Embedding", options={e.label: k for k, e in embs.items()},
                                           value=next(iter(embs)), width=210)
        self.age = pn.widgets.RadioButtonGroup(label="Age", options=["All"], value="All", color="primary",
                                               variant="outline", width=240)
        self.vmax = pn.widgets.FloatInput(label="Color max (0 = data max)", value=0.0, start=0.0, step=0.5, width=170)
        self.labels = pn.widgets.Select(label="Labels on UMAP", options=LABELS, value="", width=170)
        self.group_by = pn.widgets.Select(label="Table rows", options=GROUPS, value="cell_type", width=190)
        self.umap = UmapPlot()
        self.violin_tissue = pn.pane.Bokeh(sizing_mode="fixed")
        self.violin_cell = pn.pane.Bokeh(sizing_mode="fixed")
        # Labelled map: one fixed size for every embedding (images are scaled to fit), so it sits beside the UMAP
        # on a wide screen, wraps below otherwise, and swapping it does not re-lay out the page
        self.map = pn.pane.PNG(None, width=MAP_WIDTH, height=MAP_HEIGHT, sizing_mode="fixed")
        self.map_note = pn.pane.HTML("", width=MAP_WIDTH, sizing_mode="fixed")
        self.map_box = pn.Column(self.map_note, self.map, visible=False, width=MAP_WIDTH + 10, sizing_mode="fixed")
        umaps = pn.FlexBox(pn.pane.Bokeh(self.umap.fig, sizing_mode="fixed"), self.map_box, gap="10px",
                           sizing_mode="stretch_width")
        self.figures = pn.Tabs(("UMAP", umaps),
                               ("Top 20 tissues", self.violin_tissue), ("Top 20 cell types", self.violin_cell),
                               dynamic=True)
        self._expr: np.ndarray | None = None   # expression over the cells of the current embedding
        self.embedding.param.watch(lambda e: self._cells_changed(new_embedding=True), "value")
        self.age.param.watch(lambda e: self._cells_changed(), "value")
        self.vmax.param.watch(lambda e: self._redraw(), "value")
        self.labels.param.watch(lambda e: self.umap.set_labels(self._cells()[0], self.labels.value), "value")
        self.group_by.param.watch(lambda e: self._redraw(), "value")
        self._cells_changed(new_embedding=True)

    def panel(self):
        controls = pn.Row(self.embedding, pn.Column(pn.widgets.StaticText(value="Age"), self.age), self.vmax,
                          self.labels, sizing_mode="stretch_width")
        return pn.Column(self.header, self.card, controls, self.figures,
                         pn.Row(self.group_by, pn.Spacer(), self.download, align="end", sizing_mode="stretch_width"),
                         self.table, *self.methods, sizing_mode="stretch_width")

    # -- cells -----------------------------------------------------------------------------------------------
    def _cells(self) -> tuple[pd.DataFrame, np.ndarray]:
        cells = self.store.embeddings[self.embedding.value].cells
        if self.age.value == "All" or "age" not in cells:
            return cells, np.ones(len(cells), bool)
        mask = (cells.age == self.age.value).to_numpy()
        return cells[mask], mask

    def _cells_changed(self, new_embedding: bool = False) -> None:
        if new_embedding:
            emb = self.store.embeddings[self.embedding.value]
            self.map.object = str(emb.map) if emb.map else None
            self.map_note.object = (f"<div style='font-size:0.85em;color:#555'>Colored by "
                                    f"{esc(emb.map_colored_by.lower())}: a static picture of all cells, which does "
                                    "not follow zoom or the age filter.</div>")
            self.map_box.visible = emb.map is not None
            cells = emb.cells
            ages = [str(a) for a in cells.age.cat.categories if a != "unknown"] if "age" in cells else []
            if self.age.value not in ["All"] + ages:
                self.age.value = "All"
            self.age.options = ["All"] + ages
            self._expr = None
        cells, _ = self._cells()
        self.umap.set_cells(cells)
        self.umap.set_labels(cells, self.labels.value)
        if self.row is not None:
            self._redraw()

    # -- gene ------------------------------------------------------------------------------------------------
    def show(self, row, resolution=None):
        self.row, self.resolution = row, resolution
        self._expr = None
        self._redraw()

    def clear_plots(self) -> None:
        self._expr = None
        self.umap.clear_expression()
        self.violin_tissue.object = None
        self.violin_cell.object = None

    def _redraw(self) -> None:
        if self.row is None:
            return
        emb = self.store.embeddings[self.embedding.value]
        if self._expr is None:
            self._expr = self.store.expression(self.row, emb.key)
        cells, mask = self._cells()
        expr = self._expr[mask]
        sym = self.symbol(self.row)
        age = "" if self.age.value == "All" else f", {self.age.value}"
        self.umap.set_expression(cells, expr, self.vmax.value or None, f"{sym} expression ({emb.label}{age})")
        self.violin_tissue.object = violin_figure(expr, cells.tissue.values, title=f"{sym}: top 20 tissues by mean",
                                                  y_label="Tissue")
        self.violin_cell.object = violin_figure(expr, cells.cell_type.values,
                                                title=f"{sym}: top 20 cell types by mean",
                                                y_label="Cell ontology class")
        self._update_card(sym, emb, cells, expr)
        stats = group_stats(expr, cells[self.group_by.value])
        stats = stats.rename(columns={"group": GROUP_LABEL[self.group_by.value]})
        self.set_table(stats, f"{sym}_tabula_muris_{emb.key}_{self.group_by.value}.csv",
                       {"pct_expressing": FMT_PCT, "mean": FMT_FC, "median": FMT_FC, "q3": FMT_FC, "max": FMT_FC,
                        "cells": FMT_NUM, "expressing": FMT_NUM})

    def _update_card(self, sym: str, emb, cells: pd.DataFrame, expr: np.ndarray) -> None:
        g = self.store.genes.iloc[self.row]
        on = expr > 0
        stats = [("cells expressing", f"{int(on.sum()):,} / {len(expr):,}"),
                 ("% expressing", f"{100 * on.mean():.2f}%" if len(expr) else "NA"),
                 ("max log(CPM+1)", f"{expr.max():.2f}" if len(expr) else "NA")]
        full = self._expr
        all_cells = emb.cells
        if "age" in all_cells and self.age.value == "All":
            for a in all_cells.age.cat.categories:
                if a == "unknown":
                    continue
                m = (all_cells.age == a).to_numpy()
                stats.append((f"% expressing, {a}", f"{100 * (full[m] > 0).mean():.2f}%"))
        if on.any():
            top = group_stats(expr, cells.cell_type).iloc[0]
            stats.append(("highest mean", f"{top.group} ({top['mean']:.2f})"))
        notes = [f"Gene type: {esc(g.gene_type or 'unknown')}. Values are log(CPM + 1), natural log."]
        warnings = []
        if g.n_features > 1:
            warnings.append(f"This symbol covers {g.n_features} GENCODE vM38 genes; the quantification summed them "
                            "(transcripts were mapped to gene names).")
        if self.resolution is not None and self.resolution.note:
            notes.append(esc(self.resolution.note))
        self.card.object = card_html(sym, self.ids_text(self.row), stats, notes, warnings)
