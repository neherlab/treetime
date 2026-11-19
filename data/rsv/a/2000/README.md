# `rsv/a/2000`

## Translations

`translations/<cds>.fasta.xz` holds one amino-acid alignment per CDS of `annotation.gff3`, named by the CDS `Name` attribute. They are translations of `aln.fasta.xz` with the standard genetic code: a stop codon becomes `*`, a gap codon becomes `-`, an ambiguous or partially gapped codon becomes `X`.

```bash
./dev/docker/python python data/translate rsv/a/2000
```

## Gene alignments

`aln-G.fasta.xz` and `aln-F.fasta.xz` hold the columns of the G and F CDS of `annotation.gff3` (4652-5617 and 5697-7421) from `aln.fasta.xz`, for the example of one partition per alignment file (`timetree-gene-alignments.yaml`).

```bash
xz -dc aln.fasta.xz | awk '/^>/ {print; next} {print substr($0, 4652, 966)}' | xz -9 > aln-G.fasta.xz
xz -dc aln.fasta.xz | awk '/^>/ {print; next} {print substr($0, 5697, 1725)}' | xz -9 > aln-F.fasta.xz
```
