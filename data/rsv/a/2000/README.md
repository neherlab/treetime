# `rsv/a/2000`

## Translations

`translations/<cds>.fasta.xz` holds one amino-acid alignment per CDS of `annotation.gff3`, named by the CDS `Name` attribute. They are translations of `aln.fasta.xz` with the standard genetic code: a stop codon becomes `*`, a gap codon becomes `-`, an ambiguous or partially gapped codon becomes `X`.

```bash
./dev/docker/python python data/translate rsv/a/2000
```
