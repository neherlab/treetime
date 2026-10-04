# Taxonium JSONL output is not supported

TreeTime cannot write the JSONL format of the Taxonium tree viewer. The format holds a header line with tree metadata and a table of all mutations, then one line per node [[doc](https://github.com/theosanderson/taxonium/blob/ae4bceba28723d0b122330a7a8e1491601a97e2f/docs/taxoniumtools.md?plain=1#L39)]. Each mutation entry has a gene, a previous residue, a position, a new residue and a type, `aa` or `nt` [[src](https://github.com/theosanderson/taxonium/blob/ae4bceba28723d0b122330a7a8e1491601a97e2f/taxoniumtools/src/taxoniumtools/utils.py#L194-L216)], so the format can carry the amino-acid mutations that UShER MAT cannot ([M-io-usher-mat-rejects-amino-acid-mutations.md](M-io-usher-mat-rejects-amino-acid-mutations.md)).

Taxonium can already display TreeTime's Auspice JSON. A JSONL writer adds the format used by the UCSC daily trees and by viral_usher.

## Ecosystem

- **Producers**: `usher_to_taxonium` from `taxoniumtools` converts a MAT and a GenBank reference to JSONL. Users found by GitHub code search: the UCSC automated builds for SARS-CoV-2, mpox, RSV, dengue, influenza A and tuberculosis (`ucscGenomeBrowser/kent`, `src/hg/utils/otto/`), viral_usher (`AngieHinrichs/viral_usher`), linolium (`corbett-lab/linolium`), and Cov2Tree.org, which shows the UCSC SARS-CoV-2 tree. The search covers public default branches on GitHub only
- **Readers**: Taxonium is the only reader found
- **Maintenance**: Theo Sanderson wrote most of Taxonium and is the only author of its paper ([Sanderson 2022](https://doi.org/10.7554/eLife.82392)). Other contributors include Alex Kramer and Angie Hinrichs (UC Santa Cruz). With no schema and one main maintainer, the format can change with a `taxoniumtools` release
- **Taxodium**: the earlier name of Taxonium and of its protobuf format (`taxodium.proto` in UShER). `matUtils extract --write-taxodium` still writes it, but the current Taxonium code has no reader for it, so it is not a target

## Open questions

- The format has no schema. Which `taxoniumtools` revision defines the contract, and how is compatibility tested without running Python in the Rust tests?
- Node layout: Taxonium JSONL stores `x_dist`, optional `x_time` and `y` per node. Does TreeTime compute the layout, or must the writer match the `taxoniumtools` layout?
- Which commands write it: every command that writes MAT, or also `mugration` and `clock`?

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md): output selection design
