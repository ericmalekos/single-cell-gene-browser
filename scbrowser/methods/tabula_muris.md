#### Data
Tabula Muris and Tabula Senis FACS-sorted SmartSeq2 (plate-seq) cells from mice aged 3, 18 and 24 months, with
the consortium's cell annotations: tissue, cell ontology class, free annotation, age, sex and mouse. FASTQ files and
metadata come from the consortium's public bucket (`s3://czb-tabula-muris-senis`).

#### Quantification and normalization
- Reads were quantified with Salmon (`-l A --validateMappings --gcBias --seqBias`) against GENCODE vM38 transcripts,
  plus the spike-in and transgene sequences of the consortium's MM10-PLUS reference.
- Transcript estimates were summed per **gene name** with tximport (estimated counts, `countsFromAbundance = "no"`).
  GENCODE genes that share a name are therefore combined, and the gene card says when this applies.
- Cells with at least 500 detected genes and 50,000 counts were kept.
- Expression is log(1 + CPM), natural log, where CPM = count / cell total × 10⁶.

#### Embeddings
- **All cells, all ages** (109,429 cells): genes detected in at least one cell (75,026); 3,000 highly variable genes
  (scran `modelGeneVar`); 50 principal components; UMAP on the components (30 neighbours, min_dist 0.5; reduced
  dimension `UMAP_refined`); Louvain clustering of a shared-nearest-neighbour graph (k = 20) on the components, giving
  37 clusters. Processed object: [Zenodo 18765151](https://zenodo.org/records/18765151).
- **Immune cells, all ages** (30,411 cells of 11 cell ontology classes: B cell, CD4+ and CD8+ αβ T cell,
  granulocyte, hematopoietic stem cell, macrophage, microglial cell, monocyte, myeloid dendritic cell, neutrophil and
  NK cell): recomputed on the subset. Highly variable genes are those with a biological variance component above 0.5
  (the top 3,000 if fewer than 1,000 pass); then 50 principal components, UMAP (10 neighbours, min_dist 0.5) and
  Louvain clustering (k = 20), giving 12 clusters. Expression values are the same cells' values in the global object.
- Highly variable gene selection, PCA and UMAP were run with the random seed 42.

#### What this tab shows
- Values are the object's log(CPM + 1). A cell expresses a gene when its value is above 0.
- The table (cells, expressing, %, mean, median, upper quartile, maximum per group) and the top-20 violins (groups
  ranked by mean) are computed in the app, over the cells of the selected embedding and age.
- 661 cells have no age annotation; they appear under "All" only.
- Ensembl IDs are the GENCODE vM38 IDs of the gene name (several when genes sharing a name were combined).
- Labelled maps: the global map is drawn by the app's build from the embedding coordinates, with each tissue named at
  the median position of its cells. The immune map is the lab's figure of the same embedding.

#### References
- Tabula Muris Consortium. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. *Nature* 562,
  367–372 (2018). [doi:10.1038/s41586-018-0590-4](https://doi.org/10.1038/s41586-018-0590-4)
- Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. *Nature* 583,
  590–595 (2020). [doi:10.1038/s41586-020-2496-1](https://doi.org/10.1038/s41586-020-2496-1)
- Malekos E, Smaliy V, Carpenter S. Integrative multiomics analysis and CRISPR screening identify functional
  noncanonical translation loci in the mouse immune system. *bioRxiv* (2026).
  [doi:10.64898/2026.02.23.707583](https://doi.org/10.64898/2026.02.23.707583)
