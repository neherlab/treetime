# `sc2/2844`

2504 SARS-CoV-2 genomes sampled from 2019-12-26 to 2022-07-12, aligned to Wuhan-Hu-1 (29903 nt), from 96 countries of 6 regions.

## Files

- `tree.nwk`: the topology of the source tree of the dataset, whose origin is not recorded, with branch lengths in substitutions per site fitted by IQ-TREE 3.1.4 under GTR+G on that fixed topology, and the root next to `Wuhan-Hu-1/2019`
- `aln.fasta.xz`: one sequence per tip of the tree
- `metadata.tsv.xz`: the rows of the Nextstrain open SARS-CoV-2 metadata (`s3://nextstrain-data/files/ncov/open/metadata.tsv.zst`, version of 2026-10-18, `metadata_tsv_sha256sum` `f59776a80a1de6bf6d48acfb7f8d53005608e659098eff86c1ea09e3f2485e6d`) whose `strain`, `genbank_accession`, or `genbank_accession_rev` without its version names a tip. The `strain` column holds the tip name, so a genome that the tree names by its accession keeps that name. 4 tips have no row: `Japan/SZ-NIG-4-C4732/2022`, `OX231729`, `Switzerland/BL-ETHZ-561456/2020`, `USA/WA-S90/2020`
- `annotation.gff3`: the genes of the Wuhan-Hu-1 reference

## Regenerate

The metadata, from the open metadata of the day (rows can change between versions):

```bash
aws s3 cp s3://nextstrain-data/files/ncov/open/metadata.tsv.zst - --no-sign-request | zstd -dc > open-metadata.tsv
```

Then keep the header and every row whose `strain`, `genbank_accession`, or `genbank_accession_rev` without its version names a tip, with the tip name in `strain`.

The branch lengths, on the topology of `tree.nwk`:

```bash
xz -dkc aln.fasta.xz > aln.fasta
iqtree3 -s aln.fasta -te tree.nwk -m GTR+G -seed 42 -T 8 -o Wuhan-Hu-1/2019 --prefix fixed
```
