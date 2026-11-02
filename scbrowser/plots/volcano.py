"""Volcano plot (zeb2os_perturbseq.png style)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, LabelSet, Span
from bokeh.plotting import figure

from .common import cds_data, neglog10, style

SIG_COLOR = "#c62828"
SIG_LINE = "#5a0e0e"
NS_COLOR = "#bdbdbd"


def volcano_figure(df: pd.DataFrame, *, sig_col: str, thr: float, cutoff_p: float | None, title: str,
                   x_label: str, hover: list[tuple[str, str]], top_n: int = 15, flag_col: str | None = None,
                   highlight: str | None = None, width: int = 760, height: int = 560):
    """df needs label, logFC, PValue and sig_col; flagged points (flag_col) are drawn hollow."""
    d = df.copy()
    d["nlp"] = neglog10(d.PValue)
    d["sig"] = d[sig_col].to_numpy() < thr
    flagged = d[flag_col].astype(bool).to_numpy() if flag_col else np.zeros(len(d), bool)
    d["fill_alpha"] = np.where(flagged, 0.0, 1.0)
    # significant points filled or hollow in separate renderers: a legend key is drawn from a single data point
    ns, sg, sf = d[~d.sig], d[d.sig & ~flagged], d[d.sig & flagged]

    fig = figure(width=width, height=height, title=title, x_axis_label=x_label,
                 y_axis_label="−log10 P (edgeR QL F-test)", tools="pan,wheel_zoom,box_zoom,reset,save",
                 active_scroll="wheel_zoom", toolbar_location="above")
    r_ns = fig.scatter("logFC", "nlp", source=ColumnDataSource(cds_data(ns)), size=5, color=NS_COLOR, alpha=0.7,
                       fill_alpha="fill_alpha", line_color=NS_COLOR)
    r_sg = fig.scatter("logFC", "nlp", source=ColumnDataSource(cds_data(sg)), size=9, fill_color=SIG_COLOR,
                       line_color=SIG_LINE, line_width=0.8, legend_label=f"significant ({sig_col} < {thr:g})")
    renderers = [r_ns, r_sg]
    if len(sf):
        renderers.append(fig.scatter("logFC", "nlp", source=ColumnDataSource(cds_data(sf)), size=9, fill_alpha=0,
                                     line_color=SIG_LINE, line_width=1.2,
                                     legend_label="significant, calibration-flagged gene"))
    if highlight is not None and (d.label == highlight).any():
        own = d[d.label == highlight]
        fig.scatter("logFC", "nlp", source=ColumnDataSource(cds_data(own)), size=14, marker="diamond",
                    fill_color="#1565c0", line_color="black", legend_label=f"{highlight} itself")
    labelled = ColumnDataSource(cds_data(d[d.sig].sort_values("PValue").head(top_n)))
    fig.add_layout(LabelSet(x="logFC", y="nlp", text="label", source=labelled, x_offset=6, y_offset=2,
                            text_font_size="9pt", text_font_style="bold", text_color="#111111"))
    fig.add_tools(HoverTool(renderers=renderers, tooltips=hover))
    fig.add_layout(Span(location=0, dimension="height", line_dash="dashed", line_color="#9e9e9e", line_width=1))
    if cutoff_p:
        fig.add_layout(Span(location=float(neglog10(cutoff_p)), dimension="width", line_dash="dotted",
                            line_color="#9e9e9e", line_width=1))
    legend = fig.legend[0]
    legend.orientation = "horizontal"
    legend.label_text_font_size = "9pt"
    legend.border_line_color = None
    fig.add_layout(legend, "below")   # outside the frame, so it never covers the top hits
    if len(d):
        lim = max(0.5, float(np.nanmax(np.abs(d.logFC))) * 1.25)   # room for labels at the edges
        fig.x_range.start, fig.x_range.end = -lim, lim
        fig.y_range.start, fig.y_range.end = 0, max(1.0, float(np.nanmax(d.nlp))) * 1.08
    style(fig)
    return fig
