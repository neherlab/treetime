# `lassa/L/20`

Subsample of `lassa/L/50` (20 sequences).

> /!\ Prepared by AI.

## Methodology

1. Read metadata from source dataset
2. Parse dates from `num_date` column (numeric year or ISO format)
3. Sort sequences by date, then by `accession` for deterministic tie-breaking
4. Select 20 sequences at equal intervals across temporal range
5. Extract selected sequences from alignment using seqkit
6. Filter metadata to selected sequences
7. Build ML tree using IQ-TREE with GTR+G model, fixed seed, single thread

Date range: 1976.21 - 2022.67

## Provenance

Source checksums (MD5, first 12 chars):
- aln.fasta.xz: 0ae07e8f2e33
- metadata.tsv: f5dbd93e18bc

IQ-TREE: model=GTR+G, seed=42, threads=1

## Regenerate

See `data/subsample` script.

```bash
./dev/docker/run python3 data/subsample lassa/L/20
```

