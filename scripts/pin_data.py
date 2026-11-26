"""ONE-TIME step for the author: pin the exact input data and record it.

You do not need this to reproduce the results. It is how the files in
metadata/ were made, kept so the choice of data is transparent and can be
redone (for example with a newer ClinVar release).

What it does:
1. Chooses a dated ClinVar GRCh38 release from NCBI's archive (the latest
   in archive_2.0/<year>/, or the one given with --release YYYYMMDD) and
   reads NCBI's published md5 for it. ClinVar's "latest" file changes every
   week, so we never use it.
2. Downloads that release and finds each gene's region from ClinVar's own
   GENEINFO annotation (small variants annotated to exactly one panel gene,
   padded by --pad bases).
3. Writes metadata/regions.tsv and metadata/sources.tsv.
4. Downloads everything with scripts/download_data.py (gnomAD v4.1 joint
   slices for each gene) and writes metadata/checksums.sha256.

Then commit metadata/ together. From then on `pixi run download` only
fetches and verifies.

Usage:  pixi run pin-data                     # latest archived ClinVar
        pixi run pin-data --release 20261103   # a specific release
"""

import argparse
import csv
import datetime
import gzip
import re
import sys
import urllib.request
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import download_data as dd  # noqa: E402
from acmg_af import panel_gene  # noqa: E402

ROOT = dd.ROOT
CLINVAR_ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/{year}/"
GNOMAD_JOINT = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
                "gnomad.joint.v4.1.sites.{chrom}.vcf.bgz")
MAX_ALLELE_LENGTH = 50  # bases; larger records are structural and would stretch the regions


def fetch_text(url):
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8 (acmg-af-sarcomere)"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode()


def choose_clinvar(release):
    year = int(release[:4]) if release else datetime.date.today().year
    for y in (year, year - 1):  # early in January the new year's folder may be empty
        try:
            listing = fetch_text(CLINVAR_ARCHIVE.format(year=y))
        except OSError:
            continue  # no folder for that year (yet)
        dates = sorted(set(re.findall(r'clinvar_(\d{8})\.vcf\.gz"', listing)))
        if release:
            if release in dates:
                return CLINVAR_ARCHIVE.format(year=y) + f"clinvar_{release}.vcf.gz"
            sys.exit(f"ClinVar release {release} is not in {CLINVAR_ARCHIVE.format(year=y)}; "
                     f"available: {', '.join(dates)}")
        if dates:
            return CLINVAR_ARCHIVE.format(year=y) + f"clinvar_{dates[-1]}.vcf.gz"
    sys.exit("No ClinVar release found in the NCBI archive; give one with --release")


def provider_md5(url):
    try:
        return fetch_text(url + ".md5").split()[0]
    except OSError:
        print(f"  note: no .md5 published for {url}")
        return "-"


def gene_regions(clinvar_vcf, genes, pad):
    """Min and max position of small ClinVar variants annotated to each gene."""
    spans = {}
    with gzip.open(clinvar_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.split("\t", 8)[:8]
            if len(ref) > MAX_ALLELE_LENGTH or len(alt) > MAX_ALLELE_LENGTH:
                continue
            match = re.search(r"(?:^|;)GENEINFO=([^;]+)", info)
            gene = panel_gene(match.group(1), genes) if match else None
            if gene is None:
                continue  # not a panel gene, or two panel genes (ambiguous)
            pos = int(pos)
            lo, hi, n, _ = spans.get(gene, (pos, pos, 0, chrom))
            spans[gene] = (min(lo, pos), max(hi, pos), n + 1, chrom)
    missing = [g for g in genes if g not in spans]
    if missing:
        sys.exit(f"No ClinVar records found for {', '.join(missing)}")
    return [{"gene": g, "chrom": f"chr{spans[g][3]}", "start": max(1, spans[g][0] - pad),
             "end": spans[g][1] + pad, "clinvar_chrom": spans[g][3], "clinvar_records": spans[g][2]}
            for g in genes]


def have(source):
    """True if the file is already here and matches the provider's md5 (if any)."""
    path = ROOT / source["dest"]
    if not path.exists():
        return False
    if source["md5"] in ("", "-"):
        return path.stat().st_size > 0
    return dd.digest(path, "md5") == source["md5"]


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--release", help="ClinVar release date, YYYYMMDD (default: latest archived)")
    parser.add_argument("--pad", type=int, default=100, help="bases added on each side of a gene region")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    genes = config["genes"]

    url = choose_clinvar(args.release)
    print(f"ClinVar release: {url}")
    clinvar = [
        {"name": "clinvar_vcf", "kind": "file", "url": url, "md5": provider_md5(url),
         "dest": config["clinvar_vcf"], "index": "-"},
        {"name": "clinvar_tbi", "kind": "file", "url": url + ".tbi", "md5": provider_md5(url + ".tbi"),
         "dest": config["clinvar_vcf"] + ".tbi", "index": "-"},
    ]
    # Fetch ClinVar first: the gene regions come from it. Files already
    # downloaded and matching NCBI's md5 are kept (no second 200 MB download).
    write_tsv(dd.SOURCES, clinvar, ["name", "kind", "url", "md5", "dest", "index"])
    need = {s["dest"] for s in clinvar if not have(s)}
    dd.fetch_missing(clinvar, need)

    regions = gene_regions(ROOT / config["clinvar_vcf"], genes, args.pad)
    write_tsv(dd.REGIONS, regions, ["gene", "chrom", "start", "end", "clinvar_chrom", "clinvar_records"])
    print(f"Wrote {dd.REGIONS.relative_to(ROOT)}")

    gnomad_dir = config["gnomad_dir"]
    indexes, slices = {}, []
    for r in regions:
        remote = GNOMAD_JOINT.format(chrom=r["chrom"])
        index = f"{gnomad_dir}/index/{Path(remote).name}.tbi"
        indexes[index] = {"name": f"gnomad_index_{r['chrom']}", "kind": "file", "url": remote + ".tbi",
                          "md5": "-", "dest": index, "index": "-"}
        slices.append({"name": f"gnomad_{r['gene']}", "kind": "gnomad_slice", "url": remote,
                       "md5": "-", "dest": f"{gnomad_dir}/{r['gene']}.vcf", "index": index})
    sources = clinvar + list(indexes.values()) + slices
    write_tsv(dd.SOURCES, sources, ["name", "kind", "url", "md5", "dest", "index"])
    print(f"Wrote {dd.SOURCES.relative_to(ROOT)}")

    dd.fetch_missing(sources, {s["dest"] for s in sources if not (ROOT / s["dest"]).exists()})
    lines = [f"{dd.digest(ROOT / s['dest'], 'sha256')}  {s['dest']}" for s in sources]
    dd.CHECKSUMS.write_text("\n".join(lines) + "\n")
    dd.make_read_only([s["dest"] for s in sources])
    print(f"Wrote {dd.CHECKSUMS.relative_to(ROOT)} ({len(lines)} files). Commit metadata/ now.")


if __name__ == "__main__":
    main()
