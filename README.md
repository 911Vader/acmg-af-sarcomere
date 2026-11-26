# Silent failure in variant interpretation

How missing gnomAD frequencies distort ACMG/AMP evidence in sarcomere
cardiomyopathy genes: a reproducible analysis of the frequency criteria
(BA1, BS1, PM2) for ClinVar variants in eight sarcomere genes.

Automated variant interpretation pipelines use population allele frequencies
for three ACMG/AMP criteria. If the frequencies are never loaded (a missing
file, a wrong INFO field name, `1` vs `chr1`), the pipeline does not fail:
every variant looks absent from gnomAD, so every variant receives PM2
(evidence for pathogenicity) and BA1/BS1 can never apply. This project
measures the effect on ClinVar variants with at least two review stars in
MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1, MYL2 and MYL3: the definitive
hypertrophic cardiomyopathy genes, a core part of molecular autopsy panels
for sudden cardiac death. It also checks how robust the result is to the BS1
threshold.

It started from a real misconfiguration found in a forensic whole-genome
variant interpretation pipeline. The analysis here is a public,
self-contained reconstruction: public data only and a simplified,
independent implementation of the frequency rules, not that pipeline's code.

Author: Ghaya Alkamali (MSc Computational Biology, MBZUAI).
Capstone project for CBIO7501 Reproducible and Open Science (Fall 2026).

## Publication

None. The frequency thresholds are from Kelly et al. (2018), *Genet Med*
20:351–359, <https://doi.org/10.1038/gim.2017.218>; the criteria from
Richards et al. (2015), *Genet Med* 17:405–424,
<https://doi.org/10.1038/gim.2015.30>.

## Reproduce the results

Requirements: [Pixi](https://pixi.prefix.dev/latest/) and an internet
connection (NCBI and Google Cloud Storage must be reachable). Linux or macOS;
on Windows use WSL2 (`bcftools` has no Windows build). Nothing else to
install.

```bash
git clone https://github.com/911Vader/acmg-af-sarcomere.git
cd acmg-af-sarcomere
pixi run reproduce
```

`pixi run reproduce` installs the exact software versions in `pixi.lock`,
runs the unit tests of the rules, downloads the pinned ClinVar release and
the gnomAD slices and checks their checksums, runs the Snakemake workflow and
renders the report. A GitHub Action runs the same command on every push.

Other entry points: `pixi run test` (unit tests only, no data),
`pixi run download` (data only), `pixi run run-snakemake` (workflow),
`pixi task list` (all tasks). `pixi run pin-data` is the one-time step the
author used to choose the data and write `metadata/`. You do not need it.

To try other thresholds, genes or review-status cut-offs, edit
`config/config.yaml` and run `pixi run reproduce` again.

Exploration (not needed for the results): `pixi run explore` runs
`notebooks/01-explore-data.ipynb` headless. Optional Nextflow port of the
`classify` step (lesson 10): `pixi run -e nextflow compare-nextflow` runs it
and checks that it gives the same tables as Snakemake.

### Without Pixi

Conda or Mamba users can build the same environment from
`envs/environment.yml` (exported with `pixi workspace export
conda-environment`): `conda env create -f envs/environment.yml`, then run
`python scripts/download_data.py`, `snakemake --cores 2` and
`quarto render reports/af-acmg-report.qmd`. The versions are then not locked
as tightly as with `pixi.lock`.

### Container

The image is built from `pixi.lock`, so it has exactly the same packages:

```bash
docker build --platform linux/amd64 -f containers/Dockerfile \
  --build-arg ENV_NAME=default -t acmg-af-sarcomere:0.1.0 .
docker run --rm -v "$PWD":/work -w /work acmg-af-sarcomere:0.1.0 \
  bash -c "python scripts/download_data.py && snakemake --cores 2"
```

The *Container image* GitHub Action builds the same image on every release
and publishes it as `ghcr.io/911vader/acmg-af-sarcomere:<tag>`.
Image ID: <!-- TODO: paste the digest from the Action's summary (sha256:...) -->

## Hardware

Any laptop from the last few years. Measured on 2026-11-26 on a MacBook Pro
(Apple M5 Pro, 18 cores, 48 GB RAM, macOS 27.2): the main software
environment (`.pixi/envs/default`) takes 1.2 GB of disk, or 2.8 GB with the
optional Nextflow environment, and the raw data (`data/raw/`: the ClinVar VCF,
gnomAD indexes and slices) takes 205 MB. The Snakemake workflow ran in 21 s
from scratch, with a peak memory use under 150 MB; downloading takes a few
minutes, depending on the connection.

## Outputs

| Path | Content |
|------|---------|
| `reports/af-acmg-report.html` | The report: question, data, rules, results, threshold sensitivity, limitations, software versions |
| `results/tables/variant_evidence.tsv` | One row per ClinVar variant: gnomAD frequencies and BA1/BS1/PM2 with and without them |
| `results/tables/summary.tsv` | Counts per ClinVar label and scenario (the main result) |
| `results/tables/summary_by_gene.tsv` | The same per gene, with gnomAD frequencies |
| `results/tables/bs1_sweep.tsv` | Share of variants meeting BS1 at each threshold in `bs1_sweep` |
| `data/processed/` | The ClinVar and gnomAD tables the workflow builds |

Expected result: without frequencies, every benign variant receives PM2.
Result (ClinVar release 2026-11-15): without gnomAD, all 2,257
benign/likely benign variants and all 801 pathogenic/likely pathogenic
variants receive PM2. With gnomAD, 539 of the 2,257 benign/likely benign
variants receive BA1 (412) or BS1 (127) instead, and 0 of the 801
pathogenic/likely pathogenic variants receive benign frequency evidence
(`results/tables/summary.tsv`).

## Data

- **ClinVar**, GRCh38 VCF, a dated release from NCBI's archive
  (<https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/>).
  The exact file, NCBI's md5 and our sha256 are in `metadata/`. ClinVar data
  is freely available; see NCBI's [ClinVar terms](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
  and Landrum et al. (2018), <https://doi.org/10.1093/nar/gkx1153>.
- **gnomAD v4.1**, joint exome + genome sites VCFs
  (<https://gnomad.broadinstitute.org/downloads>), read remotely from the
  public Google Cloud bucket. Only each gene's region and five INFO fields are
  kept. gnomAD summary data are released without restriction on use
  ([terms](https://gnomad.broadinstitute.org/policies)).

Sources, gene regions and checksums are in `metadata/`. The data is not
stored in this repository. No individual-level or case data is used.

### Data availability statement

All input data are public: ClinVar (NCBI) and gnomAD v4.1 (Broad Institute).
The exact files and versions are listed with checksums in `metadata/`,
and `scripts/download_data.py` retrieves them. All code is in this
repository under the MIT license.

## Repository contents

See the table in `docs/structure.md` (each folder also has its own README).
Analysis decisions are logged in `docs/decisions.md`.

## Cite

See `CITATION.cff` (GitHub: "Cite this repository").

## License

Code: MIT (`LICENSE`). Report text and figures: CC BY 4.0
(`LICENSE-CC-BY-4.0.txt`). The input data
keep their providers' terms (see Data).

## Contact

Questions and problems: open an issue on GitHub. Contributions: see
`CONTRIBUTING.md`. Everyone taking part follows `CODE_OF_CONDUCT.md`.
