# Deposits

Put (or symlink) the deposit folders here, or point to them from `config.local.yaml`:

| folder | content | needed for |
|---|---|---|
| `erv_immdict_deposit/` | mouse Immune Dictionary counts/edgeR and human Cytokine Dictionary pseudobulks/edgeR | Mouse and Human Cytokine Dictionary tabs |
| `geiger_bmdc_deposit/` | BMDC E3-ligase Perturb-seq counts and per-target edgeR | BMDC Perturb-seq tab |

The Tabula Muris object and the annotation files are downloaded automatically by `scbrowser build`.
Everything in this folder except this README is ignored by git.
