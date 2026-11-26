# Decision log

Each entry: the decision, why, and where it is set.

- 2026-11-24: **Question and scope.** Measure what missing gnomAD frequencies
  do to the three ACMG/AMP frequency criteria only (BA1, BS1, PM2). It is not
  a full classifier: the other criteria need evidence a public dataset cannot
  give. The rules are a new, simplified implementation (`scripts/acmg_af.py`),
  not code from any production pipeline.

- 2026-11-24: **Genes.** The eight sarcomere genes with definitive evidence
  for hypertrophic cardiomyopathy (MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1,
  MYL2, MYL3). They have expert-panel frequency thresholds and are central to
  molecular autopsy after sudden cardiac death. Set in `config/config.yaml`
  (`genes`).

- 2026-11-24: **Thresholds.** The ClinGen Cardiomyopathy Expert Panel's MYH7
  values (Kelly et al. 2018): BA1 FAF ≥ 0.1 %, BS1 FAF ≥ 0.02 %, PM2 < 0.004 %,
  applied to all eight genes. A simplification: the panel later specified the
  other sarcomere genes too. The BS1 sweep in the report shows the conclusion
  does not depend on the exact value. Set in `config/config.yaml`
  (`thresholds`, `bs1_sweep`).

- 2026-11-24: **gnomAD dataset and fields.** v4.1 *joint* (exomes + genomes)
  sites, because it is what the gnomAD browser shows and covers both. BA1/BS1
  use the group-max filtering allele frequency `fafmax_faf95_max_joint`
  (ClinGen recommends the FAF); PM2 uses the group-max AF `AF_grpmax_joint`,
  falling back to `AF_joint`. Only PASS records are used
  (`gnomad_pass_only`). A missing grpmax value on a gnomAD record is read as
  "extremely rare"; gnomAD leaves it empty when allele counts are too low.

- 2026-11-24: **gnomAD slices instead of whole files.** The per-chromosome
  joint VCFs are tens of GB. `download_data.py` reads only each gene's region
  from the public bucket with bcftools and keeps five INFO fields. The slices
  are written as plain VCF with `--no-version`, so they are byte-identical on
  every run and can be checked with sha256 like any raw file. If gnomAD ever
  changes a file, the checksum check fails instead of silently giving new
  numbers.

- 2026-11-24: **ClinVar release.** A dated file from NCBI's archive, never the
  weekly "latest", which changes every week and would break the checksums.
  Chosen once by `pixi run pin-data` and recorded in `metadata/sources.tsv`.

- 2026-11-24: **Gene regions from ClinVar.** Each gene's region is the span of
  its small ClinVar variants (≤ 50 bp, annotated to exactly one panel gene)
  plus 100 bp, computed by `pin_data.py` and committed in `metadata/regions.tsv`.
  No hand-typed coordinates.

- 2026-11-24: **Reference labels.** ClinVar records with ≥ 2 stars, labelled
  Benign / Likely benign / Benign/Likely benign or Pathogenic / Likely
  pathogenic / Pathogenic/Likely pathogenic, annotated to exactly one panel
  gene. Other, non-panel genes may be listed too: ClinVar annotates every
  ACTC1 variant to the overlapping non-coding gene GJD2-DT as well (found on
  the first real run, 2026-11-26), so requiring a single gene dropped ACTC1
  entirely. Conflicting, uncertain and mixed labels (for example
  `Pathogenic|risk_factor`) are left out. Expert-panel variants are reported
  separately in the report, because the panel may have used gnomAD itself
  (circularity).

- 2026-11-24: **Matching.** ClinVar and gnomAD are joined on chromosome,
  position, REF and ALT, with no indel re-normalisation (that would need the
  3 GB reference genome). The effect is conservative: an indel represented
  differently counts as absent, which can only *add* PM2. Listed as a
  limitation in the report.

- 2026-11-25: **Rerun triggers.** The scripts read `config/config.yaml`
  directly; each Snakemake rule also lists the settings it uses as `params`.
  Changing a threshold then reruns only `classify`; changing `min_stars` reruns
  `select_clinvar` and `classify` (checked with `snakemake -n`).

- 2026-11-25: **Nextflow port (optional, lesson 10).** Only the `classify`
  step, run on the tables Snakemake built, with the same config file
  (`-params-file`). It has its own Pixi environment (`nextflow`), so the main
  environment stays small and Java is not needed to reproduce. Both
  workflows give byte-identical tables (`scripts/compare_tables.py`).

- 2026-11-25: **Container image.** Built from `pixi.lock` with the course
  Dockerfile by a GitHub Action on each release, so no one needs Docker
  installed to get it.
