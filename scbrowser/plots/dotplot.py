"""Cell type x cytokine dot plot (dotplotexample.png style).

Fill = log2 fold-change vs PBS, size = % of cells expressing the gene in that condition, thick black outline =
FDR below the threshold. A final "Baseline" group shows the PBS % expressing in gray, so its rows line up.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.models import ColorBar, ColumnDataSource, FactorRange, HoverTool, LinearColorMapper, Span
from bokeh.plotting import figure

from .common import DIVERGING, cds_data, style

COL_PX, ROW_PX = 26, 30
GROUP_PAD = 1.0
PBS_FACTOR = ("Baseline", "PBS")
UNTESTED_LINE = "#c4c4c4"
PBS_FILL = "#8c8c8c"


def dot_size(pct) -> np.ndarray:
    """Marker diameter in px; area grows with the fraction of expressing cells."""
    p = np.nan_to_num(np.asarray(pct, float), nan=0.0)
    return 4.0 + 20.0 * np.sqrt(np.clip(p, 0, 100) / 100.0)


def _separators(factors: list[tuple[str, str]]) -> list[float]:
    """Synthetic x of the gaps between top-level groups (mirrors BokehJS map_two_levels, factor_padding=0)."""
    sizes: dict[str, int] = {}
    for top, _ in factors:
        sizes[top] = sizes.get(top, 0) + 1
    out, offset = [], 0.0
    for k, n in enumerate(sizes.values()):
        offset += n
        if k < len(sizes) - 1:
            out.append(offset + GROUP_PAD / 2)
            offset += GROUP_PAD
    return out


def order_cytokines(cyts: pd.DataFrame, family_col: str) -> pd.DataFrame:
    """Keep each family contiguous: families by first appearance in the authors' order, then the order itself."""
    rank = cyts.groupby(family_col)["order"].transform("min")
    return cyts.assign(_rank=rank).sort_values(["_rank", "order"], kind="stable").drop(columns="_rank")


def dotplot_figure(long: pd.DataFrame, celltypes: list[str], cyts: pd.DataFrame, *, family_col: str,
                   fdr_thr: float, clim: float | None, title: str):
    """long: one row per (cell type, cytokine) with ct, cyt, logFC, PValue, FDR, pct, pct_pbs, n_cells, n_cells_pbs."""
    cyts = order_cytokines(cyts, family_col)
    short = family_col + "_short"
    fam_of = dict(zip(cyts.display, cyts[short]))
    factors = [(fam_of[c], c) for c in cyts.display] + [PBS_FACTOR]
    d = long[long.cyt.isin(fam_of)].copy()
    d["x"] = [(fam_of[c], c) for c in d.cyt]
    d["family"] = d.cyt.map(dict(zip(cyts.display, cyts[family_col])))
    d["size"] = dot_size(d.pct)

    tested = d[d.logFC.notna()].copy()
    if clim is None:
        top = float(np.nanmax(np.abs(tested.logFC))) if len(tested) else 1.0
        clim = max(1.0, np.ceil(top * 2) / 2)
    sig = tested.FDR.to_numpy() < fdr_thr
    tested["line_color"] = np.where(sig, "black", "#9a9a9a")
    tested["line_width"] = np.where(sig, 2.4, 0.5)
    untested = d[d.logFC.isna()]
    pbs = long.drop_duplicates("ct")[["ct", "pct_pbs", "n_cells_pbs"]].rename(
        columns={"pct_pbs": "pct", "n_cells_pbs": "n_cells"})
    pbs = pbs.assign(x=[PBS_FACTOR] * len(pbs), size=dot_size(pbs.pct))

    mapper = LinearColorMapper(palette=DIVERGING, low=-clim, high=clim)
    n_groups = len({f[0] for f in factors})
    frame = COL_PX * (len(factors) + GROUP_PAD * (n_groups - 1))
    fig = figure(x_range=FactorRange(*factors, group_padding=GROUP_PAD, factor_padding=0),
                 y_range=FactorRange(*reversed(celltypes)), width=int(frame) + 150,
                 height=ROW_PX * len(celltypes) + 260, title=title, tools="save", toolbar_location="above")

    cols = ["x", "ct", "cyt", "family", "logFC", "PValue", "FDR", "pct", "n_cells", "size"]
    src_untested = ColumnDataSource(cds_data(untested[cols]))
    src_tested = ColumnDataSource(cds_data(tested[cols + ["line_color", "line_width"]]))
    r_un = fig.scatter("x", "ct", source=src_untested, size="size", fill_color="white", line_color=UNTESTED_LINE,
                       line_width=1.0)
    r_te = fig.scatter("x", "ct", source=src_tested, size="size", fill_color={"field": "logFC", "transform": mapper},
                       line_color="line_color", line_width="line_width")
    r_pbs = fig.scatter("x", "ct", source=ColumnDataSource(cds_data(pbs)), size="size", fill_color=PBS_FILL,
                        line_color="#6e6e6e", line_width=0.5)
    fig.add_tools(
        HoverTool(renderers=[r_te], tooltips=[("Cell type", "@ct"), ("Cytokine", "@cyt"), ("Family", "@family"),
                                              ("log2FC", "@logFC{0.000}"), ("P", "@PValue{0.00e}"),
                                              ("FDR", "@FDR{0.00e}"), ("% expressing", "@pct{0.0}"),
                                              ("Cells", "@n_cells")]),
        HoverTool(renderers=[r_un], tooltips=[("Cell type", "@ct"), ("Cytokine", "@cyt"),
                                              ("Not tested", "filtered by filterByExpr"),
                                              ("% expressing", "@pct{0.0}"), ("Cells", "@n_cells")]),
        HoverTool(renderers=[r_pbs], tooltips=[("Cell type", "@ct"), ("% expressing in PBS", "@pct{0.0}"),
                                               ("Cells", "@n_cells")]),
    )
    seps = _separators(factors)
    for k, xpos in enumerate(seps):
        last = k == len(seps) - 1
        fig.add_layout(Span(location=xpos, dimension="height", line_color="black" if last else "#b0b0b0",
                            line_width=1.6 if last else 0.8))

    fig.xaxis.major_label_orientation = np.pi / 2
    fig.xaxis.group_label_orientation = "horizontal"
    fig.xaxis.group_text_font_size = "10pt"
    fig.xaxis.separator_line_alpha = 0
    fig.xgrid.grid_line_color = "#efefef"
    fig.ygrid.grid_line_color = "#efefef"
    fig.outline_line_color = None
    fig.add_layout(ColorBar(color_mapper=mapper, orientation="horizontal", title="log2 (fold-change)",
                            width=300, height=12, title_standoff=6, padding=8), "below")
    style(fig)
    fig.yaxis.major_label_text_font_size = "11pt"
    return column(fig, legend_strip(fdr_thr))


def legend_strip(fdr_thr: float):
    """Size and symbol legend (Bokeh legends draw fixed-size markers, so this is its own small figure)."""
    pcts = [1, 10, 25, 50]
    xs = list(range(len(pcts) + 3))
    fig = figure(width=110 * len(xs), height=78, x_range=(-0.5, len(xs) - 0.3), y_range=(-1.2, 1.0),
                 toolbar_location=None, tools="", min_border=0)
    fig.scatter(xs[:4], [0.2] * 4, size=dot_size(pcts), fill_color="#d9d9d9", line_color="#9a9a9a")
    fig.scatter([xs[4]], [0.2], size=14, fill_color="white", line_color="black", line_width=2.4)
    fig.scatter([xs[5]], [0.2], size=14, fill_color=PBS_FILL, line_color="#6e6e6e")
    fig.scatter([xs[6]], [0.2], size=14, fill_color="white", line_color=UNTESTED_LINE)
    labels = [f"{p}% cells" for p in pcts] + [f"FDR < {fdr_thr:g}", "PBS baseline", "not tested"]
    fig.text(xs, [-0.85] * len(xs), text=labels, text_align="center", text_font_size="9pt", text_color="#333333")
    fig.axis.visible = False
    fig.grid.visible = False
    fig.outline_line_color = None
    return fig
