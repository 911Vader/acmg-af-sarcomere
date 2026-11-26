# Silent failure in variant interpretation

Project proposal (CBIO7501 capstone). **Question:** how does a misconfigured
gnomAD annotation, where allele frequencies are never loaded, change the
ACMG/AMP frequency evidence (BA1, BS1, PM2) for ClinVar variants in eight
sarcomere cardiomyopathy genes? **Data:** a dated ClinVar GRCh38 release
(NCBI, about 200 MB, freely available) and gnomAD v4.1 joint frequencies for
the gene regions only (Broad Institute, a few MB, no restrictions on use).
**Expected output:** per-variant and summary tables of the frequency criteria
with and without gnomAD, a threshold sensitivity analysis, and a Quarto report.
