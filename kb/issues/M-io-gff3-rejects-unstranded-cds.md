# GFF3 reader rejects CDS rows without a strand

The GFF3 reader accepts only `+` and `-` in the strand column and stops on any other value [packages/treetime-io/src/gff.rs#L72-L80](../../packages/treetime-io/src/gff.rs#L72-L80). The GFF3 specification also allows `.` (not stranded) and `?` (strand unknown).

All 11 CDS rows of `data/rsv/a/20/annotation.gff3` have the strand `.`, so amino-acid reconstruction from that annotation fails before any inference:

```
treetime ancestral --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --method-anc=marginal \
  --translations='data/rsv/a/20/translations/%GENE.fasta.xz' --annotation=data/rsv/a/20/annotation.gff3 \
  --output-all=<dir>

Error:
   0: Invalid CDS strand '.' for CDS 'NS1' on GFF3 line 5
```

Listing the CDSes with `--cds-selection` instead of `--annotation` avoids the reader.

## Open question

> [!IMPORTANT]
> **Decision required.** The reader has no rule for an unstranded CDS. The GFF3 specification allows `.` and `?` in the strand column, and a bundled example annotation uses `.` on every CDS row, so `--annotation` fails on that dataset. Options:
>
> - Read `.` (and `?`) as forward strand, with or without a warning
> - Keep the error, with a message that names the accepted values and the `--cds-selection` alternative

## Smoke coverage

The `annotation` row of `dev/smoke.toml` runs the reproduction above and declares this failure.

## Related

- [M-io-gff3-cds-phase-ignored.md](M-io-gff3-cds-phase-ignored.md): another column of the same reader
