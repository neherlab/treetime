# `mpox/clade-ii/2000`

## Translations

`translations/<cds>.fasta.xz` holds one amino-acid alignment per CDS of `annotation.gff3`, named by the CDS `Name` attribute: 179 CDS of the NC_063383.1 reference, with the copies in the inverted terminal repeats named `<gene>_dup`. They are translations of `aln.fasta.xz` with the standard genetic code, on the strand of each CDS: a stop codon becomes `*`, a gap codon becomes `-`, an ambiguous or partially gapped codon becomes `X`.

```bash
./dev/docker/python python data/translate mpox/clade-ii/2000
```
