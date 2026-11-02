"""Palettes and figure styling shared by the plots."""
from __future__ import annotations

import numpy as np
from bokeh.palettes import RdBu11

FONT = "Helvetica, Arial, sans-serif"


def interp_palette(colors: list[str], n: int = 256) -> list[str]:
    """Linearly interpolate a list of hex colors to n colors."""
    rgb = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in colors], float)
    pos = np.linspace(0, 1, len(colors))
    t = np.linspace(0, 1, n)
    out = np.stack([np.interp(t, pos, rgb[:, k]) for k in range(3)], axis=1).round().astype(int)
    return ["#%02x%02x%02x" % tuple(v) for v in out]


# blue (negative) -> white -> red (positive), as in the manuscript figures; Bokeh's RdBu11 runs blue to red
DIVERGING = interp_palette(list(RdBu11))


def style(fig, title_size: str = "13pt") -> None:
    fig.title.text_font = FONT
    fig.title.text_font_size = title_size
    fig.title.text_font_style = "normal"
    for axis in fig.axis:
        axis.axis_label_text_font = FONT
        axis.axis_label_text_font_style = "normal"
        axis.axis_label_text_font_size = "11pt"
        axis.major_label_text_font = FONT
        axis.major_label_text_font_size = "10pt"
    fig.toolbar.logo = None
    fig.background_fill_color = "white"
    fig.border_fill_color = "white"


def cds_data(df) -> dict:
    """DataFrame -> ColumnDataSource data (Bokeh converts DataFrames with a reset_index pandas deprecates)."""
    return {c: (df[c].tolist() if df[c].dtype == object else df[c].to_numpy()) for c in df.columns}


def neglog10(p) -> np.ndarray:
    return -np.log10(np.clip(np.asarray(p, float), 1e-300, None))
