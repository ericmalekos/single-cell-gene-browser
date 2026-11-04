# Notes for coding agents

## Installing and building on a new machine

The steps below take a fresh clone to a running app without any files from the authors' machines. Everything is
downloaded from Zenodo and public annotation sites, and nothing needs interactive input.

Needs: Linux or macOS, conda (or mamba), internet, about 25 GB of free disk and about 16 GB of RAM. A clean-room run
from a fresh clone (Linux, 16 cores) took 8 min for the env and 1 h for `scbrowser build`, mostly downloading; it
peaked at 13.5 GB RAM (the Tabula Muris export) and left 17 GB on disk (downloads/ 6.3, data/ 3.5, app_data/ 5.0,
env 2.2).

```bash
git clone https://github.com/ericmalekos/single-cell-gene-browser.git gene_browser
cd gene_browser
conda env create -f environment.yml          # Python 3.12, R + Matrix, and the app (pip install -e .[dev])
conda activate gene_browser
scbrowser check                              # every input should read "ok" or "download"
scbrowser build                              # downloads ~10 GB (md5-checked, resumable), then builds app_data/
scbrowser check                              # every store should read "built"
python -m pytest -q                          # all pass; 1 skip (needs the authors' reference stats)
scbrowser serve --port 5006                  # http://localhost:5006/app
```

- `scbrowser build` is safe to rerun: finished downloads are kept, interrupted ones resume, and `--only <step>`
  rebuilds one step (downloads, annotations, cytokine, perturbseq, tabula_muris, tm_maps).
- With less than ~16 GB of RAM, build without Tabula Muris:
  `scbrowser build --only downloads,annotations,cytokine,perturbseq`. The other tabs work, and the Tabula Muris tab
  says "not built".
- Do not create `config.local.yaml` on a new machine unless the user already has local copies of the inputs. It is
  for pointing at existing files (see `config.local.example.yaml`).
- To check that the app serves, run `curl -s -o /dev/null -w "%{http_code}" http://localhost:5006/app`, which should
  print 200. The page itself needs a browser.
- A Mac desktop launcher is in `macos/` (see the README section "Desktop icon on a Mac").

## Repository conventions

- `config.yaml` holds the defaults and is tracked. Machine-specific paths go in `config.local.yaml`, which is
  ignored by git.
- Generated data (`app_data/`, `downloads/`, `data/*`) is never committed.
- The Methods text shown in each tab is in `scbrowser/methods/<tab>.md`.
