#### Data
The Immune Dictionary: mouse lymph-node immune cells responding to 86 cytokines or PBS in vivo (Cui et al. 2024).
Raw reads from [GEO GSE202186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202186) (44 libraries, 13
batches) were reprocessed for the ERV manuscript; cell-type labels are the authors' annotation
([Single Cell Portal SCP2554](https://singlecell.broadinstitute.org/single_cell/study/SCP2554)).

#### Reprocessing
- Reads were aligned with STARsolo (STAR 2.7.11b) to GRCm39 with GENCODE vM39 (Ensembl 116).
- Hashtags were demultiplexed with hashsolo; the cytokine of a cell comes from its hashtag call.
- Ambient RNA was removed with CellBender remove-background 0.3.2 (150 epochs, FPR 0.01) on exon-level counts, and
  doublets with scDblFinder.
- Cells kept: a single hashtag call, not a doublet, more than 1,000 UMIs, more than 500 detected genes and less than
  10% mitochondrial UMIs. This leaves 309,036 cells in the authors' 25 cell types.

#### Differential expression
- 12 cell types were tested: B cell, CD4 T, Treg, CD8 T, γδ T, NK, cDC1, cDC2, MigDC, Langerhans, pDC, macrophage.
- Per cell type: one pseudobulk of CellBender counts per mouse; samples with at least 10 cells; cytokines with at
  least 2 such mice.
- Model: design `~ batch + cytokine` with PBS as reference; `filterByExpr`; TMM normalization; `glmQLFit`; one
  quasi-likelihood F-test per cytokine against PBS (edgeR 4.8.2, R 4.5.3).
- `filterByExpr` was overridden for 9 genes, kept regardless of expression: AI506816, D17H6S56E-5 (Synir), Brip1os,
  D430020J02Rik, Gm33887, Syna, Synb, Irf7 and Stat1.
- A gene marked "not tested" in a cell type did not pass `filterByExpr` there.
- **FDR** is Benjamini-Hochberg over all genes of the fit, within each cell type and cytokine. The manuscript figures
  adjusted P values over smaller, figure-specific sets of tests, so their significance outlines can differ.

#### What this tab shows
- **Dot colour:** log2 fold change against PBS.
- **Outline:** FDR below the chosen threshold.
- **Dot size:** % of cells with at least one count, in that cell type and cytokine. The PBS column shows the PBS
  cells.
- % expressing is recomputed by the app's build from the deposited cell × gene matrix. The manuscript's dot sizes
  came from a separately filtered cell-level object; the two differ by a median of 0.07 and at most 7.3 percentage
  points.
- Cytokine names and families follow the authors' table. Mouse labels are joined to it through the deposit's
  `ligand_name_map.tsv`. Five cytokines the table does not list (AdipoQ, IFN-κ, IL-30, IL-Y and NP, all mouse only)
  are grouped under "Other".
- Ensembl IDs are the GENCODE vM39 IDs of the deposited matrix. Symbols shared by several genes get a chooser.

#### References
- Cui A, et al. Dictionary of immune responses to cytokines at single-cell resolution. *Nature* 625, 377–384 (2024).
  [doi:10.1038/s41586-023-06816-9](https://doi.org/10.1038/s41586-023-06816-9)
- Reprocessing, pseudobulks and DE: the ERV manuscript's data deposit and its README.
