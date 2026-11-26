"""Turn the gnomAD slices (one VCF per gene) into one table.

Usage: python scripts/gnomad_table.py CONFIG OUTPUT_TSV
"""

import csv
import subprocess
import sys
from pathlib import Path

import yaml


def main(config_path, output):
    config = yaml.safe_load(Path(config_path).read_text())
    fields = config["gnomad_fields"]
    names = ["faf", "af_grpmax", "af"]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t" + "\t".join(f"%INFO/{fields[n]}" for n in names) + "\n"
    rows = []
    for gene in config["genes"]:
        vcf = Path(config["gnomad_dir"]) / f"{gene}.vcf"
        out = subprocess.run(["bcftools", "query", "-f", fmt, str(vcf)],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            chrom, pos, ref, alt, filt, *values = line.split("\t")
            rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "filter": filt,
                         **dict(zip(names, values))})
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["chrom", "pos", "ref", "alt", "filter", *names],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} gnomAD records -> {output}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
