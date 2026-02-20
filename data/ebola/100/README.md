# `ebola/100`

Subsample of `ebola/362` (100 sequences).

> /!\ Prepared by AI.

## Methodology

1. Read metadata from source dataset
2. Parse dates from `" date"` column (numeric year or ISO format)
3. Sort sequences by date, then by `name` for deterministic tie-breaking
4. Select 100 sequences at equal intervals across temporal range
5. Extract selected sequences from alignment using seqkit
6. Filter metadata to selected sequences
7. Build ML tree using IQ-TREE with GTR+G model, fixed seed, single thread

Date range: 2014.21 - 2016.22

## Provenance

Source checksums (MD5, first 12 chars):
- aln.fasta.xz: 803cd57304f5
- metadata.tsv: 45c67a039d23

IQ-TREE: model=GTR+G, seed=42, threads=1

## Regenerate

See `data/subsample` script.

```bash
./dev/docker/run python3 data/subsample ebola/100
```

