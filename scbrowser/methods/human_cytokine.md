#### Data
The Parse Biosciences 10 Million Human PBMCs dataset (Oesinghaus et al. 2025). It has 9,697,974 peripheral blood
mononuclear cells from 12 donors, each treated in vitro with 90 cytokines or PBS for 24 h and profiled by
combinatorial-barcoding single-cell RNA-seq. Genes were quantified against Ensembl 109 (GENCODE v43). The data are
distributed under CC BY-NC 4.0, so reuse is restricted to non-commercial purposes.

#### Pseudobulks and differential expression
- Raw counts of the released matrix were summed per donor, cytokine and cell type, using the release's 18 cell-type
  labels and cytokine labels. Every released cell is in a pseudobulk.
- Per cell type: samples with at least 10 cells; cytokines with at least 2 donors.
- Model: design `~ donor + cytokine` with PBS as reference; `filterByExpr`; TMM normalization; `glmQLFit`; one
  quasi-likelihood F-test per cytokine against PBS (edgeR 4.8.2, R 4.5.3).
- `filterByExpr` was overridden for the 142 ERV-derived genes present in the data, kept regardless of expression.
- A gene marked "not tested" in a cell type did not pass `filterByExpr` there. ERVH48-1 (suppressyn) is not among
  the ERV-derived genes and passes only in CD4 naive T cells, so this tab shows it as tested there only. The
  deposit's repeat fits that keep it in every cell type are not loaded.
- **FDR** is Benjamini-Hochberg over all genes of the fit, within each cell type and cytokine. The manuscript figures
  adjusted P values over figure-specific sets of tests, so their FDRs differ.

#### What this tab shows
- **Dot colour:** log2 fold change against PBS.
- **Outline:** FDR below the chosen threshold.
- **Dot size:** % of cells with at least one count, in that cell type and cytokine (PBS column: the PBS cells),
  summed over donors from the deposit's per-pseudobulk counts of expressing cells.
- Over all cells, these counts equal the deposit's `human_pct_cells.tsv` exactly.
- Gene names are the release's: the Ensembl 109 name, else the Ensembl ID. Ensembl IDs come from GENCODE v43.
  Names that the release made unique with a suffix (for example MATR3-1) are separate features, offered alongside the
  base name.
- Cytokine names and families follow the authors' table.

#### References
- Oesinghaus L, et al. A single-cell cytokine dictionary of human peripheral blood. *bioRxiv* (2025).
  [doi:10.64898/2025.12.12.693897](https://doi.org/10.64898/2025.12.12.693897)
- [Parse Biosciences: 10 million human PBMCs in a single experiment](https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/)
- Pseudobulks and DE: the ERV manuscript's data deposit and its README.
