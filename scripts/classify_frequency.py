"""Apply BA1, BS1 and PM2 to each ClinVar variant, with and without gnomAD.

Two scenarios:
  with_af   frequencies from gnomAD v4.1 (joint exomes + genomes)
  no_af     the misconfiguration this project is about: the gnomAD columns
            are empty, so every variant looks "absent from controls"

Usage: python scripts/classify_frequency.py CONFIG CLINVAR_TSV GNOMAD_TSV OUTDIR

Writes, in OUTDIR:
  variant_evidence.tsv   one row per variant: frequencies, criteria in both
                         scenarios, and whether the evidence changed direction
  summary.tsv            counts per ClinVar label and scenario
  summary_by_gene.tsv    the same, per gene (with_af scenario)
  bs1_sweep.tsv          share of variants meeting BS1-or-stronger at each
                         threshold in `bs1_sweep`
"""

import math
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from acmg_af import evidence_direction, frequency_evidence, to_float  # noqa: E402

KEY = ["chrom", "pos", "ref", "alt"]


def join_gnomad(clinvar, gnomad, pass_only):
    usable = gnomad[gnomad["filter"] == "PASS"] if pass_only else gnomad
    usable = usable.drop_duplicates(KEY)
    merged = clinvar.merge(usable, on=KEY, how="left", indicator=True)
    merged["in_gnomad"] = merged.pop("_merge") == "both"
    for col in ("faf", "af_grpmax", "af"):
        merged[col] = merged[col].map(to_float)
    return merged


def apply_rules(df, thresholds):
    for scenario in ("with_af", "no_af"):
        evidence = [
            frequency_evidence(row.in_gnomad if scenario == "with_af" else False,
                               row.faf, row.af_grpmax, row.af, thresholds)
            for row in df.itertuples()
        ]
        for rule in ("BA1", "BS1", "PM2"):
            df[f"{rule}_{scenario}"] = [e[rule] for e in evidence]
        df[f"direction_{scenario}"] = [evidence_direction(e) for e in evidence]
    df["direction_changed"] = df["direction_with_af"] != df["direction_no_af"]
    return df


def summarise(df, by):
    rows = []
    for keys, part in df.groupby(by, sort=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        for scenario in ("with_af", "no_af"):
            rows.append({
                **dict(zip(by, keys)), "scenario": scenario, "variants": len(part),
                "in_gnomad": int(part["in_gnomad"].sum()) if scenario == "with_af" else 0,
                "BA1": int(part[f"BA1_{scenario}"].sum()),
                "BS1": int(part[f"BS1_{scenario}"].sum()),
                "PM2": int(part[f"PM2_{scenario}"].sum()),
                "benign_direction": int((part[f"direction_{scenario}"] == "benign").sum()),
                "pathogenic_direction": int((part[f"direction_{scenario}"] == "pathogenic").sum()),
                "no_frequency_evidence": int((part[f"direction_{scenario}"] == "none").sum()),
            })
    return pd.DataFrame(rows)


def bs1_sweep(df, values):
    rows = []
    for label, part in df.groupby("label"):
        faf = part["faf"].where(part["in_gnomad"])
        for t in values:
            meets = (faf >= t).sum()  # NaN >= t is False
            rows.append({"label": label, "threshold": t, "variants": len(part),
                         "meets_threshold": int(meets), "fraction": meets / len(part) if len(part) else math.nan})
    return pd.DataFrame(rows)


def main(config_path, clinvar_tsv, gnomad_tsv, outdir):
    config = yaml.safe_load(Path(config_path).read_text())
    clinvar = pd.read_csv(clinvar_tsv, sep="\t", dtype={"chrom": str, "ref": str, "alt": str})
    gnomad = pd.read_csv(gnomad_tsv, sep="\t", dtype=str, keep_default_na=False)
    gnomad["pos"] = gnomad["pos"].astype(int)
    # Guard against the very failure this project studies: if no gnomAD record
    # is usable, stop instead of silently treating every variant as absent.
    if config["gnomad_pass_only"] and not (gnomad["filter"] == "PASS").any():
        seen = ", ".join(sorted(gnomad["filter"].unique())[:10]) or "no records at all"
        sys.exit(f"No PASS records in the gnomAD slices (FILTER values seen: {seen}). "
                 "Check gnomad_pass_only and the slices before trusting any result.")
    df = apply_rules(join_gnomad(clinvar, gnomad, config["gnomad_pass_only"]), config["thresholds"])

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "variant_evidence.tsv", sep="\t", index=False, float_format="%.6g")
    summarise(df, ["label"]).to_csv(out / "summary.tsv", sep="\t", index=False)
    by_gene = summarise(df, ["gene", "label"])
    by_gene[by_gene["scenario"] == "with_af"].to_csv(out / "summary_by_gene.tsv", sep="\t", index=False)
    bs1_sweep(df, config["bs1_sweep"]).to_csv(out / "bs1_sweep.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"{len(df)} variants classified; {int(df['direction_changed'].sum())} change direction -> {out}/")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
