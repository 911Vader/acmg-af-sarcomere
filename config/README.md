# config

`config.yaml` holds every setting of the analysis: genes, ClinVar filters,
gnomAD fields, the BA1/BS1/PM2 thresholds and the BS1 sweep. Change a
setting here, never inside the code, then run `pixi run reproduce`.
The reason for each value is in `docs/decisions.md`.
