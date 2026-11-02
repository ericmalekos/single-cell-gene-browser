"""Labelled Tabula Muris maps: drawn from an embedding, or a configured image downscaled."""
import json

import numpy as np
import pandas as pd
from PIL import Image

from preprocess import build_tm_maps
from scbrowser import config


def test_draw_map(tmp_path):
    rng = np.random.default_rng(1)
    centers = {"Liver": (0, 0), "Lung": (8, 0), "Marrow": (0.3, 0.2), "": (4, 8)}   # "" (unannotated) is skipped
    cells = pd.concat([pd.DataFrame({"x": rng.normal(cx, 1, 200), "y": rng.normal(cy, 1, 200), "tissue": t})
                       for t, (cx, cy) in centers.items()])
    out = tmp_path / "map.png"
    build_tm_maps.draw_map(cells, "tissue", "Tissue of origin", "turbo", out)
    with Image.open(out) as im:
        assert 600 < im.width <= build_tm_maps.MAX_WIDTH and im.height > 400


def test_repel_separates_overlapping_labels():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    texts = [ax.text(5, 5, name, ha="center", va="center") for name in ("MAT", "SCAT", "GAT", "BAT")]
    fig.canvas.draw()
    build_tm_maps.repel(ax, texts)
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in texts]
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:])
    plt.close(fig)


def test_configured_image_replaces_the_drawn_map(tmp_path, monkeypatch):
    tm = tmp_path / "app_data" / "tm"
    tm.mkdir(parents=True)
    pd.DataFrame({"x": [0.0, 1, 2], "y": [0.0, 1, 0], "tissue": ["A", "B", "A"]}).to_parquet(tm / "cells_global.parquet")
    (tm / "embeddings.json").write_text(json.dumps(
        [{"key": "global", "label": "All cells", "cells": "cells_global.parquet", "expr": "expr"}]))
    Image.new("RGBA", (4000, 3000), (0, 0, 0, 0)).save(tmp_path / "mine.png")   # transparent background
    (tmp_path / "config.yaml").write_text("app_data: app_data\ntabula_muris:\n  map_images:\n    global: mine.png\n")
    monkeypatch.setenv("SCBROWSER_CONFIG", str(tmp_path / "config.yaml"))
    monkeypatch.delenv("SCBROWSER_DATA", raising=False)
    config.clear_cache()
    try:
        build_tm_maps.main()
    finally:
        config.clear_cache()
    maps = json.loads((tm / "maps.json").read_text())
    assert maps["global"]["colored_by"] == "Tissue of origin" and "mine.png" in maps["global"]["source"]
    with Image.open(tm / "map_global.png") as im:
        assert im.size == (build_tm_maps.MAX_WIDTH, build_tm_maps.MAX_WIDTH * 3 // 4)
        assert im.convert("RGB").getpixel((0, 0)) == (255, 255, 255)
