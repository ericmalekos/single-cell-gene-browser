#### Data
The E3-ligase Perturb-seq screen of Geiger-Schuller et al. in mouse bone-marrow-derived dendritic cells stimulated
with LPS for 3 h. It targets 1,130 E3 ligases, partners and substrates with 3 guides each, plus non-targeting and
cutting-control guides ([GEO GSE327057](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327057)). The
reprocessed data come from the deposit described in its README.

#### Reprocessing
- Reads were aligned with STARsolo (STAR 2.7.11b) to GRCm39 with GENCODE vM39 (Ensembl 116), 10x 3' v3 chemistry.
  Counts are GeneFull (exon + intron), with EmptyDrops_CR cell calling.
- 45 of the 46 screen channels were reprocessed, giving 504,303 cells.
- A cell was assigned to its top guide when that guide had at least 3 UMIs and at least twice the UMIs of the
  runner-up. This gives 260,195 knockout cells, 24,495 non-targeting (NTC) cells and 16,634 cutting-control cells.

#### Per-target differential expression
- For each target: one pseudobulk per channel of the target's cells, and one of that channel's NTC cells, over every
  channel that has both. GeneFull counts were summed per gene symbol.
- Model: TMM normalization; design `~ channel + condition` (NTC reference, channel-paired); `estimateDisp`,
  `glmQLFit` and `glmQLFTest` (edgeR 4.8.2, R 4.5.3). This is the model of the published WDR82 → Zeb2os analysis,
  which the reprocessing reproduces exactly.
- Every target is tested on the same 2,300 genes. These are the genes that passed `filterByExpr` on the Wdr82 vs NTC
  pseudobulks, plus 10 panel genes.
- **FDR_within_target** is Benjamini-Hochberg over the 2,300 genes of one target.
- **FDR_global** is Benjamini-Hochberg over all 1,130 × 2,300 tests. It is the one to use for "which knockouts change
  gene G".

#### Calibration
- 100 pseudo-targets, each made of 3 random non-targeting guides, were tested exactly like real targets. They contain
  no perturbation, so every call they make is false. Overall, 5.4% of their p-values are below 0.05.
- **Calibration-flagged genes** (57) reached P < 0.05 in at least 15% of the null tests. They are mostly dendritic
  cell and macrophage subpopulation markers (for example Ccl5, Cst3, Chil3 and Cd74), which a small knockout
  pseudobulk can miss, so they tend to read as "down" by chance. The KO view draws them hollow.
- The flag describes a gene, not a call: a strong effect on a flagged gene can still be real.
- The empirical FDR of the FDR_global < 0.05 calls is 9.6% over all genes and 6.7% over unflagged genes. At that
  threshold the null gives about 2 calls per target.

#### What this tab shows
- **Readout view:** one point per knockout, showing its effect on the searched gene. FDR_global is used by default,
  and the dotted line marks the largest P with FDR_global below the threshold.
- **KO view:** one point per measured gene. FDR_within_target is used by default, and the dotted line marks that
  target's largest significant P. The blue diamond is the knocked-out gene itself, when it is among the 2,300 genes.
- 59 targets are listed under outdated symbols and are also found by their current MGI symbols (for example Ddx58 is
  found by Rigi). Seven of them (for example March1, now Marchf1) are the same gene as a measured gene and are shown
  as one gene.

#### References
- Geiger-Schuller K, et al. Systematically characterizing the roles of E3-ligase family members in inflammatory
  responses with massively parallel Perturb-seq. *bioRxiv* (2023).
  [doi:10.1101/2023.01.23.525198](https://doi.org/10.1101/2023.01.23.525198)
