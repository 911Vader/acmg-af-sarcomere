# scripts

| Script | Role |
|--------|------|
| `acmg_af.py` | The frequency rules (BA1, BS1, PM2) and ClinVar helpers; pure functions, unit-tested in `tests/` |
| `download_data.py` | Downloads the raw data listed in `metadata/` and verifies checksums (`pixi run download`) |
| `pin_data.py` | One-time: chooses the ClinVar release and writes `metadata/` (`pixi run pin-data`) |
| `select_clinvar.py` | Workflow step: ClinVar records used as reference labels |
| `gnomad_table.py` | Workflow step: gnomAD slices to one table |
| `classify_frequency.py` | Workflow step: applies the rules with and without gnomAD; writes `results/tables/` |
| `compare_tables.py` | Checks that two folders of tables are identical (Snakemake vs Nextflow) |

Run scripts from the project folder through Pixi, for example
`pixi run python scripts/download_data.py --verify-only`.
