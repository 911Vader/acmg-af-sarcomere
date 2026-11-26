"""Download the raw data and check that it is exactly what we expect.

Reads   metadata/sources.tsv        where each file comes from (and the
                                    provider's md5, when the provider has one)
        metadata/regions.tsv        gene regions to cut from gnomAD
Checks  metadata/checksums.sha256   sha256 of every raw file we keep
Writes  data/raw/clinvar/           the pinned ClinVar release (VCF + index)
        data/raw/gnomad/            one gnomAD v4.1 slice per gene, plus the
                                    gnomAD indexes used to cut them
        and makes every raw file read-only (raw data is never edited).

Two kinds of source (column `kind` in sources.tsv):
  file          download the URL as is; check the provider's md5 if given
  gnomad_slice  cut the gene region (metadata/regions.tsv) out of the remote
                gnomAD VCF with bcftools, keeping only the INFO fields named
                in config/config.yaml. gnomAD's per-chromosome VCFs are tens
                of GB; the slices are a few MB.

Run it again at any time: files that are already present and correct are
not fetched again. Use --verify-only to check without downloading.
Needs bcftools on the PATH (it is in the Pixi environment).
"""

import argparse
import csv
import hashlib
import http.client
import os
import stat
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "metadata" / "sources.tsv"
REGIONS = ROOT / "metadata" / "regions.tsv"
CHECKSUMS = ROOT / "metadata" / "checksums.sha256"
CONFIG = ROOT / "config" / "config.yaml"


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def expected_sha256():
    expected = {}
    for line in CHECKSUMS.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            value, name = line.split(maxsplit=1)
            expected[name.strip().lstrip("*")] = value
    return expected


def verify(expected):
    """Return the files that are missing or have the wrong checksum."""
    bad = []
    for name, value in expected.items():
        path = ROOT / name
        if not path.exists():
            bad.append((name, "missing"))
        elif digest(path, "sha256") != value:
            bad.append((name, "checksum differs"))
    return bad


def download(url, dest, attempts=3):
    print(f"  downloading {url}", flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8 (acmg-af-sarcomere)"})
    for attempt in range(1, attempts + 1):
        try:
            # timeout: give up if the server sends nothing for 60 s
            with urllib.request.urlopen(request, timeout=60) as response, open(dest, "wb") as out:
                total = int(response.headers.get("Content-Length") or 0)
                done, shown = 0, -1
                while chunk := response.read(1 << 20):
                    out.write(chunk)
                    done += len(chunk)
                    if total > 20e6 and done * 10 // total > shown:  # progress every 10% for big files
                        shown = done * 10 // total
                        print(f"  {done / 1e6:.0f} of {total / 1e6:.0f} MB", flush=True)
            return
        except (OSError, http.client.HTTPException) as err:  # network hiccup: wait and try again
            if attempt == attempts:
                raise
            print(f"  attempt {attempt} failed ({err}); retrying in {5 * attempt} s")
            time.sleep(5 * attempt)


def run(cmd, attempts=3):
    """Run a command (bcftools reading a remote file), retrying on failure."""
    for attempt in range(1, attempts + 1):
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return
        if attempt == attempts:
            sys.exit(f"ERROR running {' '.join(cmd)}:\n{result.stderr}")
        print(f"  attempt {attempt} failed; retrying in {5 * attempt} s")
        time.sleep(5 * attempt)


def gnomad_info_fields():
    config = yaml.safe_load(CONFIG.read_text())
    fields = list(config["gnomad_fields"].values())
    return fields + [f for f in config.get("gnomad_extra_fields", []) if f not in fields]


def fetch_file(source, tmp):
    part = Path(tmp) / Path(source["dest"]).name
    download(source["url"], part)
    if source["md5"] not in ("", "-"):
        md5 = digest(part, "md5")
        if md5 != source["md5"]:
            sys.exit(f"ERROR {source['name']}: md5 {md5} does not match the provider's {source['md5']}")
        print(f"  md5 matches the provider ({md5})", flush=True)
    dest = ROOT / source["dest"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    os.replace(part, dest)


def fetch_slice(source, regions, tmp):
    """Cut one gene region out of a remote gnomAD VCF (bgzipped, tabix-indexed)."""
    gene = source["name"].removeprefix("gnomad_")
    region = regions[gene]
    index = ROOT / source["index"]
    if not index.exists():
        sys.exit(f"ERROR {source['name']}: index {source['index']} is missing (it is listed in sources.tsv)")
    keep = ",".join(f"INFO/{f}" for f in gnomad_info_fields())
    where = f"{region['chrom']}:{region['start']}-{region['end']}"
    print(f"  slicing {where} ({gene}) from {source['url']}", flush=True)
    part = Path(tmp) / f"{gene}.vcf"
    # url##idx##index: read the remote VCF with the local copy of its index.
    run(["bcftools", "view", "--no-version", "-r", where, "-O", "v", "-o", str(part) + ".full",
         f"{source['url']}##idx##{index}"])
    # Keep only the INFO fields we use: the slices stay small and readable.
    run(["bcftools", "annotate", "--no-version", "-x", f"^{keep}", "-O", "v", "-o", str(part), str(part) + ".full"])
    dest = ROOT / source["dest"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    os.replace(part, dest)


def make_read_only(names):
    for name in names:
        path = ROOT / name
        if path.exists():
            os.chmod(path, path.stat().st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def make_writable(path):
    if path.exists():
        os.chmod(path, stat.S_IWUSR | stat.S_IRUSR)


def fetch_missing(sources, names):
    """Fetch the sources whose destination is in `names` (files first, then slices)."""
    regions = {r["gene"]: r for r in read_tsv(REGIONS)} if REGIONS.exists() else {}
    (ROOT / "data").mkdir(exist_ok=True)
    # temporary files go inside the project, on the same disk as the data
    with tempfile.TemporaryDirectory(dir=ROOT / "data") as tmp:
        for kind in ("file", "gnomad_slice"):  # slices need the indexes first
            for source in sources:
                if source["kind"] != kind or source["dest"] not in names:
                    continue
                print(f"{source['name']}:", flush=True)
                make_writable(ROOT / source["dest"])
                if kind == "file":
                    fetch_file(source, tmp)
                else:
                    fetch_slice(source, regions, tmp)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verify-only", action="store_true", help="only check the files")
    args = parser.parse_args()

    if not CHECKSUMS.exists():
        sys.exit("metadata/checksums.sha256 is missing: the author must run `pixi run pin-data` once")
    expected = expected_sha256()
    bad = verify(expected)
    if bad and not args.verify_only:
        fetch_missing(read_tsv(SOURCES), {name for name, _ in bad})
        bad = verify(expected)

    for name, problem in bad:
        print(f"FAILED {name}: {problem}")
    if bad:
        sys.exit(f"{len(bad)} of {len(expected)} raw files are missing or changed")
    make_read_only(expected)
    print(f"OK: all {len(expected)} raw files present, sha256 verified, read-only")


if __name__ == "__main__":
    main()
