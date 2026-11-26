"""ACMG/AMP population-frequency rules (BA1, BS1, PM2) and ClinVar helpers.

Pure functions, no file input or output, so they can be unit-tested
(tests/test_acmg_af.py) and reused by the workflow scripts.

Only the three frequency criteria are implemented. This is not a variant
classifier: a full ACMG/AMP classification combines many more kinds of
evidence. The question here is narrower: what happens to the frequency
evidence when allele frequencies are missing?
"""

import math

# ClinVar CLNREVSTAT (VCF spelling) -> review status in gold stars
# (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,  # older releases
}


def stars(revstat):
    """Gold stars for a ClinVar CLNREVSTAT value; 0 for anything else."""
    return STARS.get(revstat, 0)


def to_float(value):
    """bcftools prints '.' for a missing value; return NaN for it."""
    if value in (None, "", "."):
        return math.nan
    return float(value)


def frequency_evidence(in_gnomad, faf, af_grpmax, af, thresholds):
    """Apply BA1, BS1 and PM2 to one variant.

    in_gnomad   the variant has a (usable) record in gnomAD
    faf         grpmax filtering allele frequency (95% CI), NaN if missing
    af_grpmax   allele frequency in the ancestry group with the highest AF
    af          overall allele frequency (used when af_grpmax is missing)
    thresholds  {"BA1": ..., "BS1": ..., "PM2": ...} as fractions

    Returns {"BA1": bool, "BS1": bool, "PM2": bool}. BS1 is reported only
    when BA1 is not met, so each variant gets at most one benign
    frequency criterion.
    """
    if not in_gnomad:
        # Absent from the population database: PM2 (supporting) applies,
        # and nothing supports a benign frequency criterion.
        return {"BA1": False, "BS1": False, "PM2": True}
    ba1 = not math.isnan(faf) and faf >= thresholds["BA1"]
    bs1 = not ba1 and not math.isnan(faf) and faf >= thresholds["BS1"]
    freq = af_grpmax if not math.isnan(af_grpmax) else af
    # gnomAD leaves grpmax empty when the allele count is too small to
    # compute it; treat that as "extremely rare" (frequency 0).
    pm2 = (0.0 if math.isnan(freq) else freq) < thresholds["PM2"]
    return {"BA1": ba1, "BS1": bs1, "PM2": pm2}


def evidence_direction(evidence):
    """Summarise the frequency evidence as one word, for tables and figures."""
    if evidence["BA1"] or evidence["BS1"]:
        return "benign"
    if evidence["PM2"]:
        return "pathogenic"
    return "none"


def label_for(clnsig, benign_labels, pathogenic_labels):
    """Map a ClinVar CLNSIG value to 'benign', 'pathogenic' or None."""
    if clnsig in benign_labels:
        return "benign"
    if clnsig in pathogenic_labels:
        return "pathogenic"
    return None


def genes_in(geneinfo):
    """Gene symbols in a ClinVar GENEINFO value, e.g. 'MYH7:4625|X:1'."""
    if geneinfo in (None, "", "."):
        return []
    return [part.split(":")[0] for part in geneinfo.split("|")]


def panel_gene(geneinfo, panel):
    """The one panel gene a ClinVar record is annotated to, else None.

    Other, non-panel genes may be listed too: every ACTC1 record is also
    annotated to the overlapping non-coding gene GJD2-DT. A record annotated
    to two panel genes is ambiguous and returns None.
    """
    hits = [g for g in genes_in(geneinfo) if g in panel]
    return hits[0] if len(hits) == 1 else None
