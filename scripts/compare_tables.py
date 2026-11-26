"""Check that two folders hold identical result tables (same names, same bytes).

Used to compare the Snakemake and Nextflow versions of the classify step.

Usage: python scripts/compare_tables.py FOLDER_A FOLDER_B
"""

import hashlib
import sys
from pathlib import Path


def digests(folder):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(folder).glob("*.tsv"))}


def main(a, b):
    da, db = digests(a), digests(b)
    if not da or not db:
        sys.exit(f"No .tsv files in {a if not da else b}: run both workflows first")
    differ = sorted(n for n in da.keys() | db.keys() if da.get(n) != db.get(n))
    for name in sorted(da.keys() | db.keys()):
        print(f"{'same ' if name not in differ else 'DIFF '} {name}")
    if differ:
        sys.exit(f"{len(differ)} table(s) differ between {a} and {b}")
    print(f"OK: all {len(da)} tables identical")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
