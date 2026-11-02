"""Figure builders, without any stores."""
import pandas as pd

from scbrowser.plots.volcano import volcano_figure


def legend_keys(fig):
    return {item.label.value: item.renderers[0].glyph for item in fig.legend[0].items}


def test_volcano_legend_keys_do_not_depend_on_the_top_point():
    df = pd.DataFrame({"label": ["Fabp4", "Rps24", "Actb", "Ccl5"], "logFC": [-2.0, -0.8, 0.1, -0.2],
                       "PValue": [1e-12, 1e-8, 0.5, 0.3], "FDR": [1e-9, 1e-6, 0.9, 0.8],
                       "flag": [True, False, False, True]})    # the most significant gene is flagged
    keys = legend_keys(volcano_figure(df, sig_col="FDR", thr=0.05, cutoff_p=None, title="", x_label="",
                                      hover=[], flag_col="flag"))
    assert keys["significant (FDR < 0.05)"].fill_alpha == 1.0
    assert keys["significant, calibration-flagged gene"].fill_alpha == 0
    # no flag column (readout view): one key, filled
    keys = legend_keys(volcano_figure(df, sig_col="FDR", thr=0.05, cutoff_p=None, title="", x_label="", hover=[]))
    assert list(keys) == ["significant (FDR < 0.05)"] and keys["significant (FDR < 0.05)"].fill_alpha == 1.0
