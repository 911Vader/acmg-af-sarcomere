// Lesson 10 (optional): the `classify` step of workflow/Snakefile, in Nextflow.
//
// It reads the same config file (-params-file config/config.yaml) and the
// tables the Snakemake workflow built, runs the same script, and writes to
// results/nextflow/tables/. `pixi run -e nextflow compare-nextflow` then checks
// that the two workflows give identical tables.
//
// Run from the project folder, after `pixi run run-snakemake`:
//   pixi run -e nextflow compare-nextflow

nextflow.enable.dsl = 2

params.clinvar_panel = "data/processed/clinvar_panel.tsv"
params.gnomad_panel = "data/processed/gnomad_panel.tsv"
params.outdir = "results/nextflow"

process CLASSIFY {
    publishDir params.outdir, mode: "copy"

    input:
    path clinvar
    path gnomad
    path config
    path scripts

    output:
    path "tables/*.tsv"

    script:
    """
    python ${scripts}/classify_frequency.py ${config} ${clinvar} ${gnomad} tables
    """
}

workflow {
    CLASSIFY(
        file(params.clinvar_panel, checkIfExists: true),
        file(params.gnomad_panel, checkIfExists: true),
        file("config/config.yaml", checkIfExists: true),
        file("scripts", type: "dir", checkIfExists: true),
    )
}
