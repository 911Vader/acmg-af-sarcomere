"""Select the ClinVar records used as reference labels.

Keeps small variants (SNVs and indels) annotated to exactly one gene of the
panel (other, non-panel genes may be listed too), with at least `min_stars` review status and a benign or pathogenic
CLNSIG (see config/config.yaml).

Usage: python scripts/select_clinvar.py CONFIG OUTPUT_TSV
"""

import csv
import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from acmg_af import label_for, panel_gene, stars  # noqa: E402

FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "INFO/GENEINFO", "INFO/CLNSIG", "INFO/CLNREVSTAT", "INFO/MC"]
COLUMNS = ["chrom", "pos", "ref", "alt", "clinvar_id", "gene", "clnsig", "label", "stars", "revstat", "consequence"]


def main(config_path, output):
    config = yaml.safe_load(Path(config_path).read_text())
    with open(config["regions"], newline="") as fh:
        regions = list(csv.DictReader(fh, delimiter="\t"))
    rows, seen = [], set()
    for region in regions:
        where = f"{region['clinvar_chrom']}:{region['start']}-{region['end']}"
        fmt = "\t".join(f"%{f}" for f in FIELDS) + "\n"
        out = subprocess.run(["bcftools", "query", "-r", where, "-f", fmt, config["clinvar_vcf"]],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            chrom, pos, vid, ref, alt, geneinfo, clnsig, revstat, mc = line.split("\t")
            if panel_gene(geneinfo, config["genes"]) != region["gene"] or alt in (".", "*") or alt.startswith("<"):
                continue
            label = label_for(clnsig, config["benign_labels"], config["pathogenic_labels"])
            n_stars = stars(revstat)
            if label is None or n_stars < config["min_stars"] or vid in seen:
                continue
            seen.add(vid)
            # MC looks like "SO:0001583|missense_variant,SO:...": keep the first term's name
            consequence = mc.split(",")[0].split("|")[-1] if mc != "." else "."
            rows.append({"chrom": f"chr{chrom}", "pos": int(pos), "ref": ref, "alt": alt, "clinvar_id": vid,
                         "gene": region["gene"], "clnsig": clnsig, "label": label, "stars": n_stars,
                         "revstat": revstat, "consequence": consequence})
    rows.sort(key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"]))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} ClinVar records selected -> {output}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
