"""Perturb-seq view: volcano of KO effects on a gene (readout) or of a KO on all genes (target), and a KO target list."""
from __future__ import annotations

import pandas as pd
import panel as pn

from ..config import dataset_info
from ..plots.volcano import volcano_figure
from ..stores import perturb_store
from .base import FMT_FC, FMT_NUM, FMT_P, DatasetView, GeneSearch, card_html, esc, fmt_p

READOUT, TARGET = "Gene as readout", "Gene as KO target"
FDR_COLUMNS = {"Global (all targets x genes)": "FDR_global", "Within target": "FDR_within_target"}
DEFAULT_FDR = {READOUT: "FDR_global", TARGET: "FDR_within_target"}
# the KO target box lists every target when clicked, so its menu needs to scroll
MENU_CSS = ".bk-menu { max-height: 340px; overflow-y: auto; }"
OVERVIEW_COLUMNS = {"target": "KO target", "current_symbol": "current symbol", "n_cells_KO": "KO cells",
                    "n_channels": "channels", "n_up_FDR_global": "genes up", "n_down_FDR_global": "genes down",
                    "n_up_FDR_within_target": "genes up", "n_down_FDR_within_target": "genes down",
                    "own_gene_logFC": "own-gene log2FC", "own_gene_FDR_within_target": "own-gene FDR",
                    "readout": "also a readout", "ti": "ti"}


class PerturbView(DatasetView):
    key = "perturbseq"

    def __init__(self):
        super().__init__(dataset_info(self.key))
        self.store = perturb_store()
        self.mode = pn.widgets.RadioButtonGroup(label="View", options=[READOUT, TARGET], value=READOUT,
                                                color="primary", variant="outline")
        self.fdr_col = pn.widgets.Select(label="FDR", options=FDR_COLUMNS, value="FDR_global", width=220)
        self.fdr = pn.widgets.Select(label="Threshold", options=[0.01, 0.05, 0.1], value=0.05, width=110)
        self.top_n = pn.widgets.Select(label="Label top", options=[0, 10, 15, 30, 50], value=15, width=100)
        self.plot = pn.pane.Bokeh(sizing_mode="fixed")
        self._roles = (-1, -1)
        self._busy = False
        self._want: str | None = None   # mode for the next gene shown (picked from the KO target list)
        self.mode.param.watch(self._mode_changed, "value")
        for w in (self.fdr_col, self.fdr, self.top_n):
            w.param.watch(lambda e: self._redraw(), "value")

        # KO target finder: a box that lists every target, and a sortable table of them with their effect counts
        self.overview = self.store.targets_overview()[list(OVERVIEW_COLUMNS)]
        named = self.overview.target + self.overview.current_symbol.map(lambda c: f" ({c})" if c else "")
        self._labels = dict(zip(self.overview.ti, named))               # target index -> "Ddx58 (Rigi)"
        self._ti_of_label = {v.lower(): k for k, v in self._labels.items()}
        self._pos_of_ti = dict(zip(self.overview.ti, range(len(self.overview))))
        n = len(self.overview)
        self.target_search = GeneSearch(label=f"KO target ({n:,} knocked out)", options=list(named), value="",
                                        case_sensitive=False, search_strategy="includes", min_characters=0,
                                        max_completions=None, restrict=True, width=300, stylesheets=[MENU_CSS],
                                        placeholder="Click for the list, or type")
        self.target_search.param.watch(self._target_typed, "value")
        flt = {"type": "input", "func": "like", "placeholder": "filter"}
        self.targets_table = pn.widgets.Tabulator(
            self.overview, titles=OVERVIEW_COLUMNS, hidden_columns=["ti"], disabled=True, selectable=1,
            show_index=False, pagination="local", page_size=10, layout="fit_data_table", sizing_mode="stretch_width",
            header_filters={"target": flt, "current_symbol": flt},
            groups={"global FDR < 0.05": ["n_up_FDR_global", "n_down_FDR_global"],
                    "within-target FDR < 0.05": ["n_up_FDR_within_target", "n_down_FDR_within_target"]},
            widths={c: 120 for c in ("n_up_FDR_global", "n_down_FDR_global", "n_up_FDR_within_target",
                                     "n_down_FDR_within_target")},
            formatters={"n_cells_KO": FMT_NUM, "own_gene_logFC": FMT_FC, "own_gene_FDR_within_target": FMT_P,
                        "readout": {"type": "tickCross", "crossElement": False}})
        self.targets_table.on_click(lambda e: self._pick_target(int(self.overview.ti.iloc[e.row])))
        both = int(((self.store.ids.readout_idx >= 0) & (self.store.ids.target_idx >= 0)).sum())
        self.finder_note = pn.pane.HTML(
            f"<div style='font-size:0.88em;color:#444;line-height:1.4'>{n:,} genes were knocked out and "
            f"{len(self.store.readouts):,} genes were measured as readouts ({both} are both). Pick a KO target "
            "here to see what its knockout changes, or search any gene in the sidebar to see which knockouts "
            "change it.</div>", sizing_mode="stretch_width", margin=(22, 10, 0, 10))

    def panel(self):
        controls = pn.Row(pn.Column(pn.widgets.StaticText(value="View"), self.mode), self.fdr_col, self.fdr,
                          self.top_n, sizing_mode="stretch_width")
        finder = pn.Row(self.target_search, self.finder_note, sizing_mode="stretch_width")
        browse = pn.Card(self.targets_table, title=f"All {len(self.overview):,} KO targets (click a row to show it)",
                         collapsed=True, sizing_mode="stretch_width")
        return pn.Column(self.header, finder, browse, self.card, controls, self.plot,
                         pn.Row(pn.Spacer(), self.download, sizing_mode="stretch_width"), self.table,
                         *self.methods, sizing_mode="stretch_width")

    def show(self, row, resolution=None):
        self.row, self.resolution = row, resolution
        self._roles = self.store.roles(row)
        gi, ti = self._roles
        available = [m for m, idx in ((READOUT, gi), (TARGET, ti)) if idx >= 0]
        want, self._want = self._want, None
        self._busy = True
        try:
            keep = next(m for m in (want, self.mode.value, available[0]) if m in available)
            self.mode.options = available
            self.mode.value = keep
            self.mode.param.trigger("value")   # new options under an unchanged value: resend the active button
            self.fdr_col.value = DEFAULT_FDR[keep]
        finally:
            self._busy = False
        self._sync_finder()
        self._redraw()

    def _sync_finder(self) -> None:
        """The KO target box and table show the target whose KO view is open, and are empty otherwise."""
        ti = self._roles[1]
        on = self.row is not None and self.mode.value == TARGET and ti >= 0
        self._busy = True
        try:
            self.target_search.value = self._labels[ti] if on else ""
            self.targets_table.selection = [self._pos_of_ti[ti]] if on else []
        finally:
            self._busy = False

    def _target_typed(self, event) -> None:
        ti = self._ti_of_label.get((event.new or "").lower())
        if not self._busy and ti is not None:
            self._pick_target(ti)

    def _pick_target(self, ti: int) -> None:
        """Show a KO target in its KO view; through the app shell, so the sidebar and other tabs follow."""
        row = int(self.store.target_rows[ti])
        if row == self.row:
            self.mode.value = TARGET
            return
        self._want = TARGET
        if self.select_gene is not None:
            self.select_gene(self.symbol(row))
        else:
            self.show(row)

    def clear_plots(self) -> None:
        self.plot.object = None
        self._sync_finder()

    def _mode_changed(self, event) -> None:
        if self._busy:
            return
        self._busy = True
        try:
            self.fdr_col.value = DEFAULT_FDR[event.new]
        finally:
            self._busy = False
        self._sync_finder()
        self._redraw()

    def _redraw(self) -> None:
        if self.row is None or self._busy:
            return
        if self.mode.value == READOUT:
            self._readout()
        else:
            self._target()

    def _readout(self) -> None:
        gi = self._roles[0]
        sym, col, thr = self.symbol(self.row), self.fdr_col.value, self.fdr.value
        df = self.store.readout_table(gi)
        sig = df[df[col] < thr]
        cutoff = self.store.cutoffs.get(thr) if col == "FDR_global" else None
        title = f"{sym} — {len(sig)} regulator{'s' if len(sig) != 1 else ''} ({df.logFC.notna().sum():,} KOs tested)"
        hover = [("KO target", "@label"), ("log2FC", "@logFC{0.000}"), ("P", "@PValue{0.00e}"),
                 ("FDR global", "@FDR_global{0.00e}"), ("FDR within target", "@FDR_within_target{0.00e}"),
                 ("KO cells", "@n_cells_KO"), ("channels", "@n_channels")]
        self.plot.object = volcano_figure(df, sig_col=col, thr=thr, cutoff_p=cutoff, title=title,
                                          x_label="log2 fold change (KO vs NT)", hover=hover, top_n=self.top_n.value)
        flag, null = bool(self.store.flag[gi]), float(self.store.null_rate[gi])
        stats = [("KO targets tested", f"{df.logFC.notna().sum():,}"),
                 (f"up, {col} < {thr:g}", str(int((sig.logFC > 0).sum()))),
                 (f"down, {col} < {thr:g}", str(int((sig.logFC < 0).sum()))),
                 ("null rate of P < 0.05", f"{100 * null:.0f}%")]
        warnings = []
        if flag:
            warnings.append(f"Calibration flag: under non-targeting pseudo-targets this gene reached P &lt; 0.05 in "
                            f"{100 * null:.0f}% of tests (≥15%). Calls on it are less reliable, though a strong "
                            "effect can still be real.")
        notes = []
        if len(sig):
            top = sig.sort_values("PValue").head(6)
            notes.append("Top regulators: " + "; ".join(
                f"{esc(r.label)} (log2FC {r.logFC:+.2f}, FDR {fmt_p(getattr(r, col))})" for r in top.itertuples()))
        emp = self.store.empirical
        if col == "FDR_global" and thr == 0.05 and len(emp):
            parts = [f"{esc(r.genes)}: {100 * r.empirical_FDR:.1f}%" for r in emp.itertuples()]
            notes.append("Empirical FDR of FDR_global &lt; 0.05 calls, from 100 non-targeting pseudo-targets: "
                         + ", ".join(parts) + ".")
        if self._roles[1] < 0:
            notes.append(f"{esc(sym)} was not knocked out in this screen, so only its readout view exists.")
        if self.resolution is not None and self.resolution.note:
            notes.append(esc(self.resolution.note))
        self.card.object = card_html(sym, self.ids_text(self.row), stats, notes, warnings)
        table = df.rename(columns={"label": "KO_target"}).sort_values("PValue", ignore_index=True)
        self.set_table(table, f"{sym}_perturbseq_readout.csv", self._formatters())

    def _target(self) -> None:
        ti = self._roles[1]
        sym, col, thr = self.symbol(self.row), self.fdr_col.value, self.fdr.value
        target = self.store.targets[ti]
        df = self.store.target_table(ti)
        sig = df[df[col] < thr]
        cutoff = self.store.cutoffs.get(thr) if col == "FDR_global" else self.store.within_target_cutoff(ti, thr)
        n_up, n_down = int((sig.logFC > 0).sum()), int((sig.logFC < 0).sum())
        title = f"{target} KO — {n_up} up, {n_down} down ({len(df):,} genes tested)"
        hover = [("Gene", "@label"), ("log2FC", "@logFC{0.000}"), ("logCPM", "@logCPM{0.0}"), ("P", "@PValue{0.00e}"),
                 ("FDR within target", "@FDR_within_target{0.00e}"), ("FDR global", "@FDR_global{0.00e}"),
                 ("calibration flag", "@calibration_flag")]
        self.plot.object = volcano_figure(df, sig_col=col, thr=thr, cutoff_p=cutoff, title=title,
                                          x_label=f"log2 fold change ({target} KO vs NT)", hover=hover,
                                          top_n=self.top_n.value, flag_col="calibration_flag", highlight=target)
        s = self.store.summary.loc[target]
        stats = [("channels", f"{int(s.n_channels)}"), ("KO cells", f"{int(s.n_cells_KO):,}"),
                 ("NTC cells", f"{int(s.n_cells_NTC):,}"),
                 (f"up / down, {col} < {thr:g}", f"{n_up} / {n_down}")]
        if pd.notna(s.own_gene_logFC):
            stats.append(("own-gene log2FC (knockdown)", f"{s.own_gene_logFC:+.2f}"))
            stats.append(("own-gene FDR within target", fmt_p(s.own_gene_FDR_within_target)))
        notes = ["Hollow points are calibration-flagged genes (null P &lt; 0.05 rate ≥ 15%), mostly subpopulation "
                 "markers that tend to read as down by chance."]
        emp = self.store.empirical
        if len(emp):
            null_calls = float(emp.iloc[0].null_calls_per_pseudotarget)
            notes.append(f"Non-targeting pseudo-targets give about {null_calls:.1f} FDR_global &lt; 0.05 calls each, "
                         "so a target with only a few calls is not distinguishable from noise.")
        if self.symbol(self.row) != target:
            notes.append(f"This KO target is listed in the screen as {esc(target)}.")
        if self._roles[0] < 0:
            notes.append(f"{esc(sym)} is not among the 2,300 readout genes, so only its KO view exists.")
        if self.resolution is not None and self.resolution.note:
            notes.append(esc(self.resolution.note))
        self.card.object = card_html(f"{sym} (KO target)", self.ids_text(self.row), stats, notes)
        table = df.rename(columns={"label": "gene"}).sort_values("PValue", ignore_index=True)
        self.set_table(table, f"{target}_KO_perturbseq.csv", self._formatters())

    @staticmethod
    def _formatters() -> dict:
        return {"logFC": FMT_FC, "logCPM": FMT_FC, "F": FMT_FC, "PValue": FMT_P, "FDR_within_target": FMT_P,
                "FDR_global": FMT_P, "n_cells_KO": FMT_NUM, "n_cells_NTC": FMT_NUM, "null_rate_P_lt_0_05": FMT_FC}
