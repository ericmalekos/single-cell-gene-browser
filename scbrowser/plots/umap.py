"""UMAP colored by expression, updated in place.

All cells are drawn once as a static background in the zero color of the colormap. Expressing cells are drawn on
top, sorted so the highest values land last, colored and enlarged by expression (as in Gaplinc_global.png). Sources
only carry numeric columns; hover labels are looked up from category lists on the client.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from bokeh.models import ColorBar, ColumnDataSource, CustomJSHover, HoverTool, LabelSet, LinearColorMapper
from bokeh.palettes import Turbo256
from bokeh.plotting import figure

from .common import style

ZERO_COLOR = Turbo256[0]
HOVER_FIELDS = [("cell_type", "Cell type"), ("tissue", "Tissue"), ("age", "Age")]
_LOOKUP_JS = "const n = names[value]; return n === undefined ? '' : n;"


class UmapPlot:
    def __init__(self, width: int = 760, height: int = 660):
        self.bg = ColumnDataSource(dict(x=np.zeros(0, np.float32), y=np.zeros(0, np.float32)))
        self.fg = ColumnDataSource(self._empty_fg())
        self.labels = ColumnDataSource(dict(x=[], y=[], text=[]))
        self.mapper = LinearColorMapper(palette=Turbo256, low=0, high=1)

        self.fig = figure(width=width, height=height, output_backend="webgl", match_aspect=True,
                          tools="pan,wheel_zoom,box_zoom,reset,save", active_scroll="wheel_zoom",
                          x_axis_label="UMAP 1", y_axis_label="UMAP 2")
        self.fig.scatter("x", "y", source=self.bg, size=2.5, color=ZERO_COLOR, alpha=0.9)
        fg = self.fig.scatter("x", "y", source=self.fg, size="size", alpha=0.9, line_color=None,
                              color={"field": "expr", "transform": self.mapper})
        self.formatters = {f"@{k}": CustomJSHover(args=dict(names=[]), code=_LOOKUP_JS) for k, _ in HOVER_FIELDS}
        tooltips = [(label, f"@{k}{{custom}}") for k, label in HOVER_FIELDS] + [("log(CPM+1)", "@expr{0.00}")]
        self.fig.add_tools(HoverTool(renderers=[fg], tooltips=tooltips, formatters=self.formatters, limit=1))
        self.color_bar = ColorBar(color_mapper=self.mapper, title="log(CPM+1)", title_standoff=8, width=14,
                                  major_label_text_font_size="10pt")
        self.fig.add_layout(self.color_bar, "right")
        self.fig.add_layout(LabelSet(x="x", y="y", text="text", source=self.labels, text_font_size="9pt",
                                     text_color="#222222", text_align="center", background_fill_color="white",
                                     background_fill_alpha=0.6))
        style(self.fig)

    @staticmethod
    def _empty_fg() -> dict:
        z = np.zeros(0, np.float32)
        return dict(x=z, y=z, expr=z, size=z, **{k: np.zeros(0, np.int32) for k, _ in HOVER_FIELDS})

    def set_cells(self, cells: pd.DataFrame) -> None:
        """New cell set (embedding or age filter): background and hover category lists.

        Axis ranges are DataRange1d with match_aspect, so they follow the background on their own."""
        self.bg.data = dict(x=cells.x.to_numpy(np.float32), y=cells.y.to_numpy(np.float32))
        for k, _ in HOVER_FIELDS:
            names = list(cells[k].cat.categories) if k in cells else []
            self.formatters[f"@{k}"].args = dict(names=[str(n) for n in names])

    def set_expression(self, cells: pd.DataFrame, expr: np.ndarray, vmax: float | None, title: str) -> None:
        on = np.flatnonzero(expr > 0)
        on = on[np.argsort(expr[on], kind="stable")]
        top = float(expr[on].max()) if len(on) else 1.0
        vmax = float(vmax) if vmax else top
        self.mapper.update(low=0, high=max(vmax, 1e-6))
        values = expr[on]
        data = dict(x=cells.x.to_numpy(np.float32)[on], y=cells.y.to_numpy(np.float32)[on],
                    expr=values.astype(np.float32),
                    size=(3.0 + 4.5 * np.clip(values / max(vmax, 1e-6), 0, 1)).astype(np.float32))
        for k, _ in HOVER_FIELDS:
            data[k] = cells[k].cat.codes.to_numpy(np.int32)[on] if k in cells else np.zeros(len(on), np.int32)
        self.fg.data = data
        self.fig.title.text = title

    def clear_expression(self) -> None:
        self.fg.data = self._empty_fg()
        self.fig.title.text = ""

    def set_labels(self, cells: pd.DataFrame, column: str | None, min_cells: int = 30) -> None:
        if not column or column not in cells:
            self.labels.data = dict(x=[], y=[], text=[])
            return
        grp = cells.groupby(column, observed=True)[["x", "y"]].median()
        counts = cells[column].value_counts()
        grp = grp[counts.reindex(grp.index).fillna(0).to_numpy() >= min_cells]
        grp = grp[grp.index.astype(str) != ""]
        self.labels.data = dict(x=grp.x.to_numpy(), y=grp.y.to_numpy(), text=[str(i) for i in grp.index])
