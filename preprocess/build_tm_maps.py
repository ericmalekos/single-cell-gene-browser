"""Labelled reference maps for the Tabula Muris tab: one static PNG per embedding.

    python -m preprocess.build_tm_maps      (or `scbrowser build --only tm_maps`)

Writes app_data/tm/map_<embedding>.png and maps.json. The global embedding is colored by tissue of origin and the
immune subset by cell type, each group named at its median position, in the style of the lab's figures (turbo for
tissues, magma for cell types). An image set under tabula_muris.map_images in the config is used instead
(downscaled); it must be drawn from the same embedding, or it will not line up with the expression UMAP.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from scbrowser.config import app_data, load_config, resolve

# embedding -> (cell column, legend title, colormap)
STYLE = {"global": ("tissue", "Tissue of origin", "turbo"), "immune": ("cell_type", "Cell type", "magma")}
MAX_WIDTH = 1400   # pixels: twice the ~700 the app shows, for sharp text on high-density screens


def log(msg: str) -> None:
    print(msg, flush=True)


def repel(ax, texts: list, pad: float = 3.0, iterations: int = 300) -> None:
    """Move overlapping labels apart (a small ggrepel, in display pixels); long moves get a leader line."""
    renderer = ax.figure.canvas.get_renderer()
    to_px, to_data = ax.transData, ax.transData.inverted()
    anchor = np.array([to_px.transform(t.get_position()) for t in texts])
    pos = anchor.copy()
    size = np.array([[b.width, b.height] for b in (t.get_window_extent(renderer) for t in texts)])
    for _ in range(iterations):
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                d = pos[j] - pos[i]
                over = (size[i] + size[j]) / 2 + pad - np.abs(d)
                if (over > 0).all():                 # boxes overlap: split along the axis that needs less
                    k = 1 if over[1] <= over[0] else 0
                    step = (np.sign(d[k]) or 1.0) * over[k] / 2
                    pos[i, k] -= step
                    pos[j, k] += step
                    moved = True
        if not moved:
            break
    for t, p, a, (_, h) in zip(texts, pos, anchor, size):
        t.set_position(to_data.transform(p))
        if np.hypot(*(p - a)) > 1.5 * h:
            (ax0, ay0), (ax1, ay1) = to_data.transform(a), to_data.transform(p)
            ax.plot([ax0, ax1], [ay0, ay1], color="#333333", lw=0.6, zorder=2)


def draw_map(cells: pd.DataFrame, column: str, title: str, cmap: str, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patheffects import withStroke

    labels = cells[column].astype(str).to_numpy()
    keep = labels != ""
    x, y, labels = cells.x.to_numpy()[keep], cells.y.to_numpy()[keep], labels[keep]
    groups = sorted(set(labels), key=str.lower)
    colors = plt.get_cmap(cmap)(np.linspace(0, 1, len(groups)))
    codes = pd.Categorical(labels, categories=groups).codes
    order = np.random.default_rng(0).permutation(len(x))   # no group hidden under another by drawing order

    fig, ax = plt.subplots(figsize=(6.2, 5.6), dpi=220)
    ax.scatter(x[order], y[order], c=colors[codes[order]], s=1.2 if len(x) < 50_000 else 0.5, linewidths=0,
               rasterized=True)
    ax.set_aspect("equal")
    ax.margins(0.03)
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    med = pd.DataFrame({"x": x, "y": y, "g": labels}).groupby("g")[["x", "y"]].median()
    texts = [ax.text(r.x, r.y, g, fontsize=8.5, fontweight="bold", ha="center", va="center", zorder=3,
                     path_effects=[withStroke(linewidth=2.2, foreground="white")]) for g, r in med.iterrows()]
    handles = [Line2D([], [], marker="o", ls="", color=c, markersize=6) for c in colors]
    ax.legend(handles, groups, title=title, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              ncol=2 if len(groups) > 14 else 1, fontsize=8, title_fontsize=9, alignment="left",
              handletextpad=0.2, columnspacing=0.8, labelspacing=0.3)
    fig.canvas.draw()   # fixes the axis limits and text sizes that repel() works with
    repel(ax, texts)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    shrink(out, out)


def shrink(src: Path, out: Path, max_width: int = MAX_WIDTH) -> None:
    """Copy an image as a PNG no wider than max_width, on white, with a 256-color palette.

    The app embeds the maps in the page; the palette makes them 3-5 times smaller (they have few distinct colors)."""
    from PIL import Image
    with Image.open(src) as im:
        rgba = im.convert("RGBA")
    flat = Image.new("RGB", rgba.size, "white")
    flat.paste(rgba, mask=rgba.getchannel("A"))   # transparent figure backgrounds become white
    if flat.width > max_width:
        flat = flat.resize((max_width, round(flat.height * max_width / flat.width)), Image.LANCZOS)
    small = flat.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE)
    palette = small.getpalette()
    for i in range(0, len(palette), 3):
        if min(palette[i:i + 3]) >= 250:          # the octree leaves white at 254
            palette[i:i + 3] = [255, 255, 255]
    small.putpalette(palette)
    small.save(out, optimize=True)


def main() -> None:
    root = app_data() / "tm"
    embeddings = json.loads((root / "embeddings.json").read_text())
    images = load_config()["tabula_muris"].get("map_images") or {}
    maps = {}
    for e in embeddings:
        key = e["key"]
        column, title, cmap = STYLE.get(key, ("cell_type", "Cell type", "turbo"))
        out = root / f"map_{key}.png"
        image = resolve(images[key]) if images.get(key) else None
        if image is not None:
            shrink(image, out)
            source = f"image {image.name}"
        else:
            draw_map(pd.read_parquet(root / e["cells"], columns=["x", "y", column]), column, title, cmap, out)
            source = "drawn from the embedding"
        maps[key] = dict(file=out.name, colored_by=title, source=source)
        log(f"{out.name}: {title.lower()}, {source} ({out.stat().st_size / 1e3:.0f} KB)")
    (root / "maps.json").write_text(json.dumps(maps, indent=2))


if __name__ == "__main__":
    main()
