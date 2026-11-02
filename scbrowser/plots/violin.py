"""Top-N group violins with jittered points (cell_violin.png / tissue_violin.png style)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, FactorRange, HoverTool
from bokeh.plotting import figure

from .common import style

HALF_WIDTH = 0.42


def _density(values: np.ndarray, grid: np.ndarray, bw: float) -> np.ndarray:
    # binned Gaussian KDE: exact enough for display and independent of the group size
    step = grid[1] - grid[0]
    hist, _ = np.histogram(values, bins=len(grid), range=(grid[0] - step / 2, grid[-1] + step / 2))
    k = np.arange(-int(4 * bw / step) - 1, int(4 * bw / step) + 2) * step
    kernel = np.exp(-0.5 * (k / bw) ** 2)
    return np.convolve(hist, kernel / kernel.sum(), mode="same")


def violin_figure(values: np.ndarray, groups: pd.Categorical, *, title: str, top: int = 20,
                  x_label: str = "log(CPM + 1)", y_label: str = "", max_points: int = 300,
                  zero_points: int = 40, seed: int = 0, width: int = 720):
    df = pd.DataFrame({"g": groups, "v": values})
    df = df[df.g.astype(str) != ""]
    stats = df.groupby("g", observed=True).v.agg(n="size", mean="mean", n_expr=lambda s: int((s > 0).sum()))
    stats = stats.sort_values(["mean", "n_expr"], ascending=False).head(top)
    order = [str(g) for g in stats.index]
    xmax = max(float(df.v.max()) if len(df) else 1.0, 1.0)
    grid = np.linspace(0, xmax * 1.02, 256)
    bw = max(0.08 * xmax, 0.1)

    patches = dict(xs=[], ys=[], group=[], n=[], pct=[], mean=[])
    points = dict(x=[], y=[])
    rng = np.random.default_rng(seed)
    labels, vals = df.g.astype(str).to_numpy(), df.v.to_numpy()
    for g in order:
        v = vals[labels == g]
        dens = _density(v, grid, bw)
        dens = dens / dens.max() * HALF_WIDTH if dens.max() > 0 else dens
        keep = dens > 1e-3 * HALF_WIDTH
        gx, gd = grid[keep], dens[keep]
        patches["xs"].append(np.concatenate([gx, gx[::-1]]).tolist())
        patches["ys"].append([(g, float(d)) for d in gd] + [(g, -float(d)) for d in gd[::-1]])
        patches["group"].append(g)
        patches["n"].append(len(v))
        patches["pct"].append(100 * float((v > 0).mean()))
        patches["mean"].append(float(v.mean()))
        nz = v[v > 0]
        if len(nz) > max_points:
            nz = rng.choice(nz, max_points, replace=False)
        zeros = np.zeros(min(zero_points, int((v == 0).sum())))
        pts = np.concatenate([zeros, nz])
        points["x"].extend(pts.tolist())
        points["y"].extend((g, float(o)) for o in rng.uniform(-0.12, 0.12, len(pts)))

    fig = figure(width=width, height=max(260, 34 * len(order) + 90), title=title, x_axis_label=x_label,
                 y_axis_label=y_label, y_range=FactorRange(*reversed(order)), tools="save,reset",
                 toolbar_location="above")
    pr = fig.patches("xs", "ys", source=ColumnDataSource(patches), fill_color="#d0d0d0", line_color="#b5b5b5",
                     fill_alpha=0.9)
    fig.scatter("x", "y", source=ColumnDataSource(points), size=3.2, color="black", alpha=0.6)
    fig.add_tools(HoverTool(renderers=[pr], tooltips=[("Group", "@group"), ("Cells", "@n"),
                                                     ("% expressing", "@pct{0.0}"), ("Mean", "@mean{0.000}")]))
    fig.x_range.start = -0.02 * xmax
    fig.ygrid.grid_line_color = "#ececec"
    style(fig)
    return fig
