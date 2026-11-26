# Project structure

| Folder | What goes in it | Rule |
|--------|-----------------|------|
| `data/raw/` | ClinVar release and gnomAD slices, downloaded by `scripts/download_data.py` | Read-only, not committed |
| `data/processed/` | ClinVar and gnomAD tables built by the workflow | Rebuilt, not committed |
| `metadata/` | Data sources, gene regions, checksums | Committed |
| `config/` | All analysis settings (`config.yaml`) | Committed |
| `scripts/` | Rules module, download and pinning scripts, workflow steps | Committed |
| `tests/` | Unit tests of the frequency rules | Committed |
| `workflow/` | `Snakefile` (the workflow); `main.nf` (optional Nextflow port of `classify`) | Committed |
| `reports/` | Quarto report source and bibliography | Source committed, HTML not |
| `results/` | Tables and logs made by the workflow | Rebuilt, not committed |
| `notebooks/` | Exploration notebook (`pixi run explore`) | Committed, outputs cleared |
| `docs/` | This file, the decision log, software versions, workflow graph | Committed |
| `envs/` | Conda export of the Pixi environment | Generated |
| `containers/` | `Dockerfile` and `apptainer.def`, built from `pixi.lock` | Committed |
| `.github/workflows/` | `reproduce.yml` (clean-clone check), `container.yml` (image on release), `pin-data.yml` (one-time data pinning) | Committed |

How the pieces connect:

```
metadata/sources.tsv ─┐
metadata/regions.tsv ─┼─ scripts/download_data.py ─► data/raw/{clinvar,gnomad}/
metadata/checksums ───┘                                   │
config/config.yaml ──────────────────────────────────────►│
                                                          ▼
workflow/Snakefile:  select_clinvar ─┐
                     gnomad_table ───┴─► classify ─► results/tables/*.tsv
                                                          │
reports/af-acmg-report.qmd ◄──────────────────────────────┘ ─► reports/af-acmg-report.html
```
