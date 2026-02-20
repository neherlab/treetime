# `zika/20`

Subsample of `zika/86` (20 sequences).

> /!\ Prepared by AI.

## Methodology

1. Read metadata from source dataset
2. Parse dates from `date` column (numeric year or ISO format)
3. Sort sequences by date, then by `#name` for deterministic tie-breaking
4. Select 20 sequences at equal intervals across temporal range
5. Extract selected sequences from alignment using seqkit
6. Filter metadata to selected sequences
7. Build ML tree using IQ-TREE with GTR+G model, fixed seed, single thread

Date range: 2013.82 - 2016.45

## Provenance

Source checksums (MD5, first 12 chars):
- aln.fasta.xz: b51b4ab87d0a
- metadata.tsv: adeace34ba1d

IQ-TREE: model=GTR+G, seed=42, threads=1

## Regenerate

See `data/subsample` script.

```bash
./dev/docker/run python3 data/subsample zika/20
```

