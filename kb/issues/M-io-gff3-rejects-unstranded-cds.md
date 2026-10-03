# GFF3 reader rejects CDS rows without a strand

The GFF3 reader accepts only `+` and `-` in the strand column and stops on any other value [packages/treetime-io/src/gff.rs#L100-L108](../../packages/treetime-io/src/gff.rs#L100-L108). The GFF3 specification also allows `.` (not stranded) and `?` (strand unknown).

Every CDS row of `data/rsv/a/20/annotation.gff3` has the strand `.`, so amino-acid reconstruction from that annotation fails before any inference:

```
treetime ancestral --tree=data/rsv/a/20/tree.nwk --aln=data/rsv/a/20/aln.fasta.xz --method-anc=marginal \
  --translations='data/rsv/a/20/translations/%GENE.fasta.xz' --annotation=data/rsv/a/20/annotation.gff3 \
  --output-all=<dir>

Error:
   0: Invalid CDS strand '.' for CDS 'NS1' on GFF3 line 5
```

Listing the CDSes with `--cdses` instead of `--annotation` avoids the reader.

## Open question

How should an unstranded CDS be read?

- As forward strand, with or without a warning
- As an error, as now, with a message that names the accepted values

## Smoke coverage

The `annotation` row of `dev/smoke.toml` runs the reproduction above and declares this failure.

## Related

- [M-io-gff3-cds-phase-ignored.md](M-io-gff3-cds-phase-ignored.md): another column of the same reader
