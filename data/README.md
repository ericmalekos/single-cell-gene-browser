# Deposits

`scbrowser build` downloads the deposit files it needs into these folders from Zenodo. To use copies you already
have, put or symlink them here, or point to them from `config.local.yaml`.

| folder | Zenodo record | content | needed for |
|---|---|---|---|
| `erv_immdict_deposit/` | [23896882](https://doi.org/10.5281/zenodo.23896882) | mouse Immune Dictionary counts/edgeR and human Cytokine Dictionary pseudobulks/edgeR | Mouse and Human Cytokine Dictionary tabs |
| `geiger_bmdc_deposit/` | [23896808](https://doi.org/10.5281/zenodo.23896808) | BMDC E3-ligase Perturb-seq counts and per-target edgeR | BMDC Perturb-seq tab |

The Tabula Muris objects and the annotation files are downloaded to `downloads/`.
Everything in this folder except this README is ignored by git.
