# workflow

`Snakefile`: the analysis steps (`select_clinvar`, `gnomad_table`,
`classify`), run with `pixi run run-snakemake`. Settings come from
`config/config.yaml`. The graph of the steps is in `docs/structure.md`;
`pixi run rulegraph` writes Snakemake's own version to `docs/rulegraph.dot`.

`main.nf`: optional Nextflow port of the `classify` step (lesson 10), run with
`pixi run -e nextflow compare-nextflow` after the Snakemake workflow.
