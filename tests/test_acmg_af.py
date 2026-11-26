"""Unit tests of the frequency rules. Run with: pixi run test"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from acmg_af import evidence_direction, frequency_evidence, genes_in, label_for, panel_gene, stars, to_float  # noqa: E402

T = {"BA1": 0.001, "BS1": 0.0002, "PM2": 0.00004}
NAN = math.nan


def test_absent_from_gnomad_gets_pm2_only():
    assert frequency_evidence(False, NAN, NAN, NAN, T) == {"BA1": False, "BS1": False, "PM2": True}


def test_missing_frequencies_look_like_absence():
    # The misconfiguration: values exist in reality but were never loaded.
    common = frequency_evidence(False, 0.05, 0.08, 0.06, T)  # in_gnomad=False hides them
    assert common == {"BA1": False, "BS1": False, "PM2": True}


def test_common_variant_gets_ba1_not_bs1():
    assert frequency_evidence(True, 0.002, 0.003, 0.001, T) == {"BA1": True, "BS1": False, "PM2": False}


def test_thresholds_are_inclusive():
    assert frequency_evidence(True, 0.001, 0.002, 0.001, T)["BA1"]
    assert frequency_evidence(True, 0.0002, 0.0003, 0.0001, T)["BS1"]


def test_between_bs1_and_ba1():
    assert frequency_evidence(True, 0.0005, 0.0007, 0.0004, T) == {"BA1": False, "BS1": True, "PM2": False}


def test_rare_variant_in_gnomad_still_gets_pm2():
    assert frequency_evidence(True, 0.000001, 0.00002, 0.00001, T) == {"BA1": False, "BS1": False, "PM2": True}


def test_pm2_falls_back_to_overall_af_when_grpmax_missing():
    assert not frequency_evidence(True, NAN, NAN, 0.0001, T)["PM2"]
    assert frequency_evidence(True, NAN, NAN, NAN, T)["PM2"]  # no counts at all: extremely rare


def test_neither_benign_nor_pm2():
    # Above the PM2 cut-off but below BS1: no frequency criterion applies.
    e = frequency_evidence(True, 0.00005, 0.0001, 0.00005, T)
    assert e == {"BA1": False, "BS1": False, "PM2": False}
    assert evidence_direction(e) == "none"


def test_direction():
    assert evidence_direction({"BA1": True, "BS1": False, "PM2": False}) == "benign"
    assert evidence_direction({"BA1": False, "BS1": True, "PM2": False}) == "benign"
    assert evidence_direction({"BA1": False, "BS1": False, "PM2": True}) == "pathogenic"


def test_stars():
    assert stars("reviewed_by_expert_panel") == 3
    assert stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert stars("criteria_provided,_single_submitter") == 1
    assert stars("no_assertion_criteria_provided") == 0


def test_labels():
    benign, pathogenic = ["Benign", "Likely_benign"], ["Pathogenic"]
    assert label_for("Likely_benign", benign, pathogenic) == "benign"
    assert label_for("Pathogenic", benign, pathogenic) == "pathogenic"
    assert label_for("Uncertain_significance", benign, pathogenic) is None


def test_parsing():
    assert genes_in("MYH7:4625") == ["MYH7"]
    assert genes_in("MYH7:4625|MYH6:4624") == ["MYH7", "MYH6"]
    assert genes_in(".") == []
    assert math.isnan(to_float("."))
    assert to_float("0.5") == 0.5


def test_panel_gene():
    panel = ["MYH7", "ACTC1", "TPM1"]
    assert panel_gene("MYH7:4625", panel) == "MYH7"
    # ClinVar annotates every ACTC1 variant to the overlapping GJD2-DT as well
    assert panel_gene("ACTC1:70|GJD2-DT:101928174", panel) == "ACTC1"
    assert panel_gene("ACTC1:70|TPM1:7168", panel) is None  # two panel genes: ambiguous
    assert panel_gene("BRCA2:675", panel) is None
    assert panel_gene(".", panel) is None
