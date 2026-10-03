# Taxonium JSONL output is not supported

TreeTime cannot write the JSONL format of the Taxonium tree viewer. The format holds a header line with tree metadata and a table of all mutations, then one line per node [[doc](https://github.com/theosanderson/taxonium/blob/ae4bceba28723d0b122330a7a8e1491601a97e2f/docs/taxoniumtools.md?plain=1#L39)]. Each mutation entry has a gene, a previous residue, a position, a new residue and a type, `aa` or `nt` [[src](https://github.com/theosanderson/taxonium/blob/ae4bceba28723d0b122330a7a8e1491601a97e2f/taxoniumtools/src/taxoniumtools/utils.py#L194-L216)], so the format can carry the amino-acid mutations that UShER MAT cannot ([M-io-usher-mat-rejects-amino-acid-mutations.md](M-io-usher-mat-rejects-amino-acid-mutations.md)).

Taxonium can already display TreeTime's Auspice JSON. A JSONL writer adds the format used by the UCSC daily trees and by viral_usher.

## Open questions

- The format has no schema. Which `taxoniumtools` revision defines the contract, and how is compatibility tested without running Python in the Rust tests?
- Node layout: Taxonium JSONL stores `x_dist`, optional `x_time` and `y` per node. Does TreeTime compute the layout, or must the writer match the `taxoniumtools` layout?
- Which commands write it: every command that writes MAT, or also `mugration` and `clock`?

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md): output selection design
