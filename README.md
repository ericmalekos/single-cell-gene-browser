# Single-cell gene browser

A web app for looking up one gene across four single-cell datasets. Type a gene symbol, alias or Ensembl ID
(ENSMUSG…/ENSG…, with or without version). Each dataset tab then shows that gene's summary, a figure, and a
sortable table you can download as CSV.

| Tab | Data | Figure |
|---|---|---|
| Tabula Muris | Tabula Muris / Senis FACS SmartSeq2, requantified against GENCODE vM38 ([Zenodo 18765151](https://zenodo.org/records/18765151)): 109,429 cells, plus an immune subset of 30,411 cells ([Zenodo 23926318](https://doi.org/10.5281/zenodo.23926318)) | UMAP colored by log(CPM+1), beside a labelled map (tissue of origin, or cell type for the immune subset); top-20 tissue / cell-type violins |
| Mouse Immune Dictionary | Immune Dictionary (GSE202186) reprocessed; edgeR vs PBS for 12 cell types ([Zenodo 23896882](https://doi.org/10.5281/zenodo.23896882)) | cell type × cytokine dot plot |
| Human Cytokine Dictionary | Parse 10M PBMC; edgeR vs PBS for 18 cell types (same Zenodo record) | the same dot plot |
| BMDC Perturb-seq | E3-ligase Perturb-seq (GSE327057) reprocessed; per-target edgeR ([Zenodo 23896808](https://doi.org/10.5281/zenodo.23896808)) | volcano, with the gene as readout or as KO target; a list of the 1,130 KO targets |

## Install

You need conda (or Python ≥ 3.11).

```bash
git clone https://github.com/ericmalekos/single-cell-gene-browser.git gene_browser
cd gene_browser
conda env create -f environment.yml      # installs the app, the `scbrowser` command and R
conda activate gene_browser
```

Without conda, use a virtual environment instead (building the Tabula Muris tab then needs R with the `Matrix`
package on the PATH):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Get the data

The app reads one folder of prebuilt stores, `app_data/` (about 5 GB). There are two ways to get it.

### A. Use a prebuilt copy

This option needs no deposits and no R. Copy an `app_data/` folder from someone who has built it, then point to it
in one of three ways:

- put the folder in the repository root, or
- `export SCBROWSER_DATA=/path/to/app_data`, or
- set `app_data: /path/to/app_data` in `config.local.yaml`.

Then skip to [Run](#run).

### B. Build from the Zenodo deposits

1. **Inputs.** `scbrowser build` downloads the data below, checking each file's md5 (R comes with the conda environment):

   | input | Zenodo record | needed for |
   |---|---|---|
   | Immune / Cytokine Dictionary deposit (1.8 GB of it) | [23896882](https://doi.org/10.5281/zenodo.23896882) | the two cytokine tabs |
   | BMDC Perturb-seq deposit (1.9 GB of it) | [23896808](https://doi.org/10.5281/zenodo.23896808) | the Perturb-seq tab |
   | Tabula Muris object (5.1 GB) | [18765151](https://zenodo.org/records/18765151) | the Tabula Muris tab |
   | Tabula Muris immune subset (1.3 GB) | [23926318](https://doi.org/10.5281/zenodo.23926318) | its second UMAP embedding |
   | GENCODE vM38 and v43 GTFs, MGI and HGNC symbol tables | (from GENCODE, MGI and HGNC) | Ensembl ID and alias search |
   | R with the `Matrix` package | installed by `environment.yml` (or system R) | Tabula Muris export only |

   The deposit files go to `data/`; the others go to `downloads/`. The human files of the cytokine deposit are
   derived from the Parse dataset and carry its CC BY-NC 4.0 non-commercial term.

2. **Optional: use copies you already have.** Copy `config.local.example.yaml` to `config.local.yaml` and set
   their paths there. The build then downloads only what is missing. In `config.local.yaml` you can also:
   - point to existing copies of the Tabula Muris objects (to skip those downloads);
   - leave out the immune subset (`tabula_muris.immune_rds: null`);
   - use your own labelled UMAP figures instead of the drawn maps (`tabula_muris.map_images`). Each figure must
     show the same embedding as the app, or it will not line up with the expression UMAP.

3. **Check, then build.**

   ```bash
   scbrowser check      # shows which inputs were found and which stores exist
   scbrowser build      # downloads what is missing, then builds every store it has inputs for
   ```

   Once the inputs are downloaded, a complete build takes about 5 minutes. Expect:
   - about 20 GB of disk: the downloads (about 10 GB) and `app_data/`;
   - about 25 GB of RAM during the Tabula Muris export.

   Steps whose inputs are missing are skipped, and the app shows those tabs as "not built". To rebuild one step,
   use `scbrowser build --only cytokine` (the steps are downloads, annotations, cytokine, perturbseq,
   tabula_muris and tm_maps). `tm_maps` takes a few seconds: it draws the labelled maps from the Tabula Muris
   store, global cells by tissue of origin and immune cells by cell type.

   Every step checks itself against the deposits:
   - DE values equal the deposit panel tables;
   - human expressing-cell counts equal `human_pct_cells.tsv` exactly;
   - mouse % expressing is within the deposit README's stated deviation;
   - Wdr82 → Zeb2os is logFC +1.818;
   - the immune subset's values are identical to the global object's.

## Run

```bash
scbrowser serve --show                    # http://localhost:5006/app
```

To let others on the network use it, run `scbrowser serve --address 0.0.0.0`. It prints the URLs. If people reach
the machine by a name it doesn't know, add `--origin <that-name>`. Any other option is passed to `panel serve`, for
example `--num-procs 4`.

Links can carry the gene and the tab: `http://localhost:5006/app?gene=Synir&tab=1` (tabs are numbered from 0 in the
table above).

## Desktop icon on a Mac

The lab Mac needs only the prebuilt data. It does not need the deposits or R.

1. **Copy the folder.** Copy this folder, including `app_data/` (about 5 GB), to the Mac. You can use an external
   drive, or run this from the Linux machine:
   `rsync -a --exclude downloads/ --exclude config.local.yaml gene_browser/ <user>@<mac>:gene_browser/`
   For rsync, first turn on Remote Login on the Mac (System Settings → General → Sharing).
2. **Run the installer.** On the Mac, run `bash ~/gene_browser/macos/install_mac.sh` in Terminal. It:
   - installs Miniforge if conda is missing;
   - creates the environment, which takes a few minutes;
   - adds **Gene Browser** to `~/Applications`, with a shortcut on the Desktop.

**Using it:** double-click the icon and the browser opens at http://localhost:5006/app. The first start takes up to
about 30 seconds. Double-click again to reopen the page or to stop the server. The log is
`~/Library/Logs/GeneBrowser.log`.

- **Data elsewhere** (for example on an external drive): `bash macos/install_mac.sh --data /Volumes/Drive/app_data`.
- **Update:** copy the new code (and `app_data/`, if it was rebuilt), then rerun the installer.
- **Remove:** delete the app and the Desktop shortcut, then run `conda env remove -n gene_browser`.
- **Change the icon:** replace `macos/icon.png` with any 1024 × 1024 PNG, then rerun the installer.

## Tests

```bash
pytest
```

Tests that need built stores are skipped until they exist.

## How a query is resolved

The lookup tries these in order, within the active dataset:

1. Ensembl ID, with versions and `_PAR_Y` ignored.
2. Exact symbol.
3. Case-insensitive symbol.
4. Alias.
5. Otherwise, suggestions.

Aliases come from MGI and HGNC synonyms plus `aliases.tsv` (for example `Synir` → D17H6S56E-5 and `suppressyn` →
ERVH48-1); add your own rows there.

The query is kept when you switch tabs, so `Stat1` finds `STAT1` in the human data. That is a match by name, not an
orthology mapping, and the app says so. An Ensembl ID from the other species is translated to its symbol first.

Where several genes match, a chooser appears in the sidebar. This happens with duplicated mouse symbols and with
anndata `-1` variants such as `MATR3-1` in the human data.

Several Perturb-seq targets use outdated symbols, so current symbols find them: `Rigi` finds Ddx58 and `Marchf1`
finds March1. A target whose current symbol is also a readout gene is merged with that readout.

On the Perturb-seq tab, the **KO target** box lists only the 1,130 knocked-out genes, with current symbols in
brackets (click it for the whole list, or type part of a name). Below it, a table of all KO targets can be sorted,
for example by how many genes each knockout changes. Picking a target from either one opens its KO view and
searches for it in every tab.

## What the numbers are

Each tab ends with a collapsed **Methods** section, written from the deposit READMEs and the Tabula Muris pipeline
(README and scripts). The text is in `scbrowser/methods/<tab>.md`; edit it there.

- **Cytokine FDR** is the deposit's genome-wide Benjamini-Hochberg adjustment (all genes within each cell type ×
  cytokine). The manuscript figures adjusted over smaller test sets, so their outlines can differ.
- **% expressing:**
  - Mouse values are recomputed for all genes from the deposited cell × gene matrix.
  - The manuscript tables cover 14 genes and come from a slightly different cell set.
  - Human values are summed `n_expressing` over donors.
- **Perturb-seq significance:**
  - The readout view uses `FDR_global`, as the deposit README recommends for "which perturbations change gene G".
  - The KO view uses `FDR_within_target`.
  - Hollow points are calibration-flagged genes.
- **Tabula Muris expression** is log(CPM + 1), natural log.
  - Symbols that cover several GENCODE genes were summed during quantification, and the gene card says so.
  - 661 cells have no age annotation. They appear under "All" only.
- **Licence:** the human data derive from a CC BY-NC 4.0 dataset, so reuse is restricted to non-commercial purposes.

## Troubleshooting

- **`scbrowser build` skips Tabula Muris.** R or its Matrix package is missing; install them as shown in the inputs
  table above.
- **The page stays blank when opened from another machine.** The browser's host name is not an allowed websocket
  origin. Restart with `--origin <host name>`.
- **A tab says "not built".** Run `scbrowser check` to see which input is missing, then `scbrowser build --only <step>`.

## Layout

```
config.yaml                defaults (tracked); config.local.yaml holds machine-specific paths (ignored)
data/                      deposits, downloaded by `scbrowser build` or placed by you (ignored)
macos/                     Mac installer, launcher and icon
preprocess/                builders behind `scbrowser build`, the R export and resources/
scbrowser/                 the app: cli, config, gene index, stores, plots/, views/, methods/ (tab texts), ui.py, app.py
tests/                     pytest suite
app_data/, downloads/      generated (ignored)
```
