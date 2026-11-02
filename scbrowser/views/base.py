"""Shared parts of a dataset view: header, gene card, results table with CSV download."""
from __future__ import annotations

import html
import io
import re
from pathlib import Path
from typing import Callable

import pandas as pd
import panel as pn
import param
from bokeh.models.widgets.tables import NumberFormatter, ScientificFormatter

from ..genes import GeneIndex, Resolution

FMT_FC = NumberFormatter(format="0.000")
FMT_PCT = NumberFormatter(format="0.0")
FMT_NUM = NumberFormatter(format="0,0")
# p-values and FDRs always in scientific notation (the default switches to fixed below 1e-3 and prints 0)
FMT_P = ScientificFormatter(precision=2, power_limit_low=0, power_limit_high=1)
METHODS_DIR = Path(__file__).resolve().parent.parent / "methods"

CARD_CSS = """
.gb-card { font-family: Helvetica, Arial, sans-serif; color: var(--panel-on-background-color, #222); }
.gb-card h3 { margin: 0 0 2px 0; font-size: 1.35em; }
.gb-card .ids { color: #666; font-size: 0.9em; margin-bottom: 8px; word-break: break-all; }
.gb-card .stats { display: flex; flex-wrap: wrap; gap: 8px; margin: 6px 0; }
.gb-card .stat { border: 1px solid #ddd; border-radius: 6px; padding: 6px 10px; min-width: 120px; background: #fafafa; }
.gb-card .stat .v { font-size: 1.15em; font-weight: 600; }
.gb-card .stat .k { font-size: 0.8em; color: #666; }
.gb-card .note { font-size: 0.86em; color: #555; margin: 3px 0; }
.gb-card .warn { font-size: 0.88em; color: #8a5300; background: #fff4e0; border-radius: 4px; padding: 4px 8px;
                 margin: 4px 0; }
.gb-header { font-size: 0.9em; color: #444; line-height: 1.35; }
.gb-header .cite { color: #777; font-size: 0.92em; }
.gb-header .licence { color: #8a5300; font-size: 0.92em; }
"""


class GeneSearch(pn.widgets.AutocompleteInput):
    """AutocompleteInput with Bokeh's max_completions, which Panel does not expose."""
    max_completions = param.Integer(default=20, allow_None=True)


def esc(x) -> str:
    return html.escape(str(x))


def card_html(title: str, ids: str = "", stats: list[tuple[str, str]] = (), notes: list[str] = (),
              warnings: list[str] = ()) -> str:
    parts = [f"<style>{CARD_CSS}</style><div class='gb-card'><h3>{esc(title)}</h3>"]
    if ids:
        parts.append(f"<div class='ids'>{esc(ids)}</div>")
    parts += [f"<div class='warn'>{w}</div>" for w in warnings]
    if stats:
        parts.append("<div class='stats'>" + "".join(
            f"<div class='stat'><div class='v'>{esc(v)}</div><div class='k'>{esc(k)}</div></div>" for k, v in stats)
            + "</div>")
    parts += [f"<div class='note'>{n}</div>" for n in notes]
    parts.append("</div>")
    return "".join(parts)


def methods_card(key: str) -> list[pn.layout.Card]:
    """Collapsed card with scbrowser/methods/<key>.md (none if missing); links open in a new tab, keeping the app."""
    path = METHODS_DIR / f"{key}.md"
    if not path.exists():
        return []
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2" target="_blank" rel="noopener">\1</a>',
                  path.read_text(encoding="utf-8"))
    return [pn.Card(pn.pane.Markdown(text, sizing_mode="stretch_width", max_width=1100), title="Methods",
                    collapsed=True, sizing_mode="stretch_width")]


def fmt_p(p: float) -> str:
    return "NA" if pd.isna(p) else (f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}")


class DatasetView:
    key = ""
    store = None

    def __init__(self, info: dict):
        self.info = info
        lines = [f"<div class='gb-header'>{esc(info.get('description', ''))}",
                 f"<div class='cite'>{esc(info.get('citation', ''))}</div>"]
        if info.get("licence"):
            lines.append(f"<div class='licence'>{esc(info['licence'])}</div>")
        self.header = pn.pane.HTML(f"<style>{CARD_CSS}</style>" + "".join(lines) + "</div>",
                                   sizing_mode="stretch_width")
        self.card = pn.pane.HTML("", sizing_mode="stretch_width")
        self.methods = methods_card(self.key)   # views add it at the bottom of their panel
        self.table = pn.widgets.Tabulator(pd.DataFrame(), pagination="local", page_size=12, disabled=True,
                                          show_index=False, sizing_mode="stretch_width", layout="fit_data_table")
        self.download = pn.widgets.FileDownload(callback=self._csv, filename="table.csv", label="Download table (CSV)",
                                                color="primary", variant="outline", width=230)
        self._table_df = pd.DataFrame()
        self.row: int | None = None
        self.resolution: Resolution | None = None
        self.select_gene: Callable[[str], None] | None = None   # set by the app shell: show a gene in every tab

    @property
    def index(self) -> GeneIndex:
        return self.store.index

    def symbol(self, row: int) -> str:
        return self.index.symbols[row]

    def ids_text(self, row: int) -> str:
        ids = self.index.ids[row]
        return ", ".join(ids) if ids else "no Ensembl ID"

    def set_table(self, df: pd.DataFrame, filename: str, formatters: dict) -> None:
        self._table_df = df
        self.table.formatters = formatters
        self.table.value = df
        self.download.filename = filename

    def _csv(self) -> io.BytesIO:
        # utf-8-sig so spreadsheet programs keep the Greek cytokine names
        return io.BytesIO(self._table_df.to_csv(index=False).encode("utf-8-sig"))

    def clear(self, message: str) -> None:
        self.row = None
        self.card.object = card_html("No gene selected", notes=[esc(message)])
        self.set_table(pd.DataFrame(), "table.csv", {})
        self.clear_plots()

    def clear_plots(self) -> None:
        """Remove the previous gene's figures (views override)."""

    def show(self, row: int, resolution: Resolution | None = None) -> None:
        raise NotImplementedError

    def panel(self) -> pn.viewable.Viewable:
        raise NotImplementedError
