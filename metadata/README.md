# metadata

Files that describe the raw data. They are small and committed. They were
written once by `pixi run pin-data` (`scripts/pin_data.py`); edit them only
by running that script again.

## sources.tsv

Where each raw file comes from. Columns:

| Column | Meaning |
|--------|---------|
| `name` | short name of the source |
| `kind` | `file` (download as is) or `gnomad_slice` (cut a gene region out of a remote gnomAD VCF with bcftools) |
| `url` | download link, or the remote VCF for a slice |
| `md5` | checksum published by the provider (NCBI's `.md5` files); `-` where the provider publishes none |
| `dest` | where the file goes in the project |
| `index` | for slices: the local copy of the remote VCF's tabix index |

## regions.tsv

One row per gene: `chrom`, `start`, `end` (GRCh38, 1-based, inclusive; the
span of the gene's small ClinVar variants plus 100 bp), `clinvar_chrom` (the
same chromosome as ClinVar names it, without `chr`) and `clinvar_records`
(small ClinVar records that defined the span).

## checksums.sha256

SHA-256 of every raw file we keep, including the gnomAD slices. Check with
`sha256sum -c metadata/checksums.sha256`, or `pixi run download`, which
checks before and after downloading.
