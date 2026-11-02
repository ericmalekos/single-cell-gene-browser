"""App shell: shared gene search in the sidebar, one tab per dataset, views built on first use."""
from __future__ import annotations

import traceback

import panel as pn

from .config import dataset_info
from .genes import Resolution
from .views.base import GeneSearch, card_html, esc

DATASETS = ["tabula_muris", "mouse_cytokine", "human_cytokine", "perturbseq"]
DEFAULT_GENE = "Irf7"


def make_view(key: str):
    if key == "tabula_muris":
        from .views.tabula_muris import TabulaMurisView
        return TabulaMurisView()
    if key in ("mouse_cytokine", "human_cytokine"):
        from .views.cytokine import CytokineView
        return CytokineView(key.split("_")[0])
    from .views.perturbseq import PerturbView
    return PerturbView()


class BrowserApp:
    def __init__(self):
        self.search = GeneSearch(label="Gene", value=DEFAULT_GENE, case_sensitive=False, search_strategy="starts_with",
                                 min_characters=2, restrict=False, max_completions=20, sizing_mode="stretch_width",
                                 placeholder="Symbol or Ensembl ID, then Enter")
        self.status = pn.pane.HTML("", sizing_mode="stretch_width")
        self.chooser = pn.widgets.Select(label="Several genes match; showing", visible=False,
                                         sizing_mode="stretch_width")
        self.suggestions = pn.FlexBox(sizing_mode="stretch_width")
        self.placeholders = {k: pn.Column(sizing_mode="stretch_width", min_height=400) for k in DATASETS}
        self.tabs = pn.Tabs(*[(dataset_info(k)["title"], self.placeholders[k]) for k in DATASETS],
                            sizing_mode="stretch_width")
        self.views: dict = {}
        self.stale = set(DATASETS)
        self._choosing = False

        self.search.param.watch(lambda e: self._query_changed(), "value")
        self.tabs.param.watch(lambda e: self._activate(), "active")
        self.chooser.param.watch(self._chose, "value")

        help_text = pn.pane.HTML(
            "<div style='font-size:0.85em;color:#555;line-height:1.4'>Search by symbol (any case), alias or "
            "Ensembl ID, with or without version. The query carries across tabs; a mouse symbol matches its "
            "human counterpart by name, which is not an orthology mapping.</div>", sizing_mode="stretch_width")
        # no dark-theme toggle: the figures and cards are drawn for a white background
        self.template = pn.template.FastListTemplate(
            title="Single-cell gene browser", sidebar_width=330, accent="#2b5f9e", theme_toggle=False,
            sidebar=[self.search, self.status, self.chooser, self.suggestions, pn.layout.Divider(), help_text],
            main=[self.tabs])
        if pn.state.location is not None:
            pn.state.location.sync(self.search, {"value": "gene"})
            pn.state.location.sync(self.tabs, {"active": "tab"})
        pn.state.onload(self._activate)

    # -- events ----------------------------------------------------------------------------------------------
    @property
    def active_key(self) -> str:
        return DATASETS[self.tabs.active]

    def _view(self, key: str):
        if key not in self.views:
            holder = self.placeholders[key]
            holder.loading = True
            try:
                view = make_view(key)
                view.select_gene = self.select_gene
                self.views[key] = view
                holder.objects = [view.panel()]
            except FileNotFoundError as exc:   # store not built: say how to build it instead of failing the page
                title = esc(dataset_info(key)["title"])
                holder.objects = [pn.pane.Alert(
                    f"<b>{title} has not been built.</b> Run <code>scbrowser build</code> (see the README for the "
                    f"inputs it needs), then reload this page.<br><small>Missing: {esc(exc.filename or exc)}</small>",
                    alert_type="warning")]
                return None
            except Exception as exc:
                holder.objects = [pn.pane.Alert(f"Could not load {key}: {esc(exc)}<pre>{esc(traceback.format_exc())}"
                                                "</pre>", alert_type="danger")]
                return None
            finally:
                holder.loading = False
        return self.views[key]

    def _activate(self) -> None:
        key = self.active_key
        view = self._view(key)
        if view is None:
            return
        self.search.options = view.index.completions()
        if key in self.stale:
            self._apply(view)

    def select_gene(self, symbol: str) -> None:
        """Search for a gene picked inside a view (the Perturb-seq KO target list)."""
        if self.search.value == symbol:
            self._query_changed()
        else:
            self.search.value = symbol

    def _query_changed(self) -> None:
        self.stale = set(DATASETS)
        view = self.views.get(self.active_key)
        if view is not None:
            self._apply(view)

    def _apply(self, view) -> None:
        self.stale.discard(view.key)
        res = view.index.resolve(self.search.value)
        if not res.found:
            self._not_found(view, res)
            return
        self._choosing = True
        try:
            self.chooser.options = {view.index.label(r): r for r in res.rows}
            self.chooser.value = res.rows[0]
            self.chooser.visible = len(res.rows) > 1
        finally:
            self._choosing = False
        self.suggestions.objects = []
        self.status.object = self._status_html(view, res, res.rows[0])
        self._show(view, res.rows[0], res)

    def _chose(self, event) -> None:
        if self._choosing or event.new is None:
            return
        view = self.views.get(self.active_key)
        if view is not None:
            res = view.resolution
            self.status.object = self._status_html(view, res, event.new)
            self._show(view, event.new, res)

    def _show(self, view, row: int, res: Resolution | None) -> None:
        try:
            view.show(row, res)
        except Exception as exc:
            view.card.object = card_html("Error", warnings=[esc(f"{type(exc).__name__}: {exc}")])
            traceback.print_exc()

    def _not_found(self, view, res: Resolution) -> None:
        self.chooser.visible = False
        q = esc(res.query) if res.query else "nothing"
        title = esc(view.info["title"])
        reason = f" ({esc(res.how)})" if res.how else ""
        self.status.object = (f"<div style='color:#8a5300'><b>{q}</b> was not found in {title}{reason}.</div>"
                              + ("<div style='font-size:0.85em;color:#555'>Did you mean:</div>"
                                 if res.suggestions else ""))
        buttons = []
        for s in res.suggestions:
            b = pn.widgets.Button(label=s, color="light", height=30, margin=(2, 3))
            b.on_click(lambda e, s=s: setattr(self.search, "value", s))
            buttons.append(b)
        self.suggestions.objects = buttons
        view.clear(f"{res.query!r} is not in this dataset." if res.query else "Type a gene symbol or Ensembl ID.")

    @staticmethod
    def _status_html(view, res: Resolution | None, row: int) -> str:
        how = f" <span style='color:#777'>by {esc(res.how)}</span>" if res and res.how else ""
        note = f"<div style='font-size:0.85em;color:#8a5300'>{esc(res.note)}</div>" if res and res.note else ""
        return f"<div><b>{esc(view.index.label(row))}</b>{how}</div>{note}"


def build_app():
    return BrowserApp().template
