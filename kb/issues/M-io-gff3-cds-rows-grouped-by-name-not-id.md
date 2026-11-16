# GFF3 CDS rows are grouped by name, not by `ID`

The GFF3 reader joins `CDS` rows into one CDS by the first name attribute it finds, `Name` before `ID` (`NAME_ATTRS_CDS` and `fn resolve_cds_name()` in [packages/treetime-io/src/gff.rs](../../packages/treetime-io/src/gff.rs)), ignores `Parent`, and sorts the segments of each CDS by coordinate ([packages/treetime-io/src/gff.rs#L93](../../packages/treetime-io/src/gff.rs#L93)).

The GFF3 specification joins the rows of one feature by a shared `ID`, and `Name` need not be unique ([GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md?plain=1)). Two consequences follow:

- Two CDS with the same `Name` and different `ID` merge into one CDS whose segments are the rows of both
- The relation between a `gene` and its `CDS` rows is lost, which [M-io-gff3-gene-only-annotations-yield-no-cds.md](M-io-gff3-gene-only-annotations-yield-no-cds.md) needs

Nextclade groups adjacent rows that share an `ID`, builds the hierarchy from `Parent`, keeps segments in file order, and wraps circular genomes ([gene.rs#L42-L130](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/gene.rs#L42), [cds.rs#L46-L90](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/cds.rs#L46-L90)).

> [!IMPORTANT]
> **Investigation required.** Sorting by coordinate gives the biological order of the segments only for forward-strand CDS without a wrap. Check how the readers of the parsed segments order minus-strand CDS with several segments, and CDS across the origin of a circular genome, before relying on the sort.

## Fix direction

Group rows by `ID` and `Parent`, with the name as a fallback for rows without `ID`; keep the segments in biological 5' to 3' order; use the `Name` only as the CDS name; stop with an error when two CDS share a name, because the outputs key CDS by name. Axis 6 of [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) describes the whole reader.

## Validation

- Two CDS with one `Name` and different `ID` give an error that names both
- One CDS with two rows that share an `ID` gives one CDS with two segments in file order
- The annotations of `data/` give the same CDS as before
