# GFF3 annotations with only `gene` rows yield no CDS

The GFF3 reader skips every row whose type is not `CDS` ([packages/treetime-io/src/gff.rs#L56-L58](../../packages/treetime-io/src/gff.rs#L56-L58)). An annotation that describes its coding regions with `gene` rows only gives no CDS.

- `data/sc2/2844/annotation.gff3` and `data/sc2/5100/annotation.gff3` have 12 `gene` rows and no `CDS` row
- 38 of 120 virus annotations (117 in the Nextclade datasets repository plus the 3 distinct annotations in `data/`) have `gene` rows only, among them SARS-CoV-2 `orfs`, dengue, yellow fever, and most influenza HA and NA files

Nextclade, which produces the annotations of the Nextstrain datasets, uses a `gene` without `CDS` children as a CDS ([gene.rs#L42-L130](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/gene.rs#L42)).

## Fix direction

Read `gene` rows, build the feature hierarchy from `ID` and `Parent` ([M-io-gff3-cds-rows-grouped-by-name-not-id.md](M-io-gff3-cds-rows-grouped-by-name-not-id.md)), and use a `gene` without `CDS` children as a CDS. Report each such fallback in the run record ([M-cli-run-warnings-and-diagnostics-only-on-stderr.md](M-cli-run-warnings-and-diagnostics-only-on-stderr.md)). Axis 6 of [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) describes the whole reader.

## Validation

- `--annotation data/sc2/5100/annotation.gff3` gives one CDS per `gene` row, with the coordinates and strand of the row
- A `gene` with `CDS` children gives the CDS of its children only
