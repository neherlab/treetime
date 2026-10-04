# The multi-format tree I/O decision describes removed converters

[kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) describes a `convert` subcommand, a `ConverterGraph` intermediate representation, and graph converters for Auspice, UShER MAT and PhyloXML (`AuspiceRead`, `UsherRead`, `PhyloxmlToGraph` and their write counterparts). Its links point to `packages/treetime-cli/src/convert/` and `packages/treetime-io/src/phyloxml.rs`, which no longer exist.

These parts were removed:

- `ebfcc094` refactor(app-cli): remove convert binary
- `b944399c` refactor: remove unreachable usher-read, initialize_partitions, reseed_transitional_from_payloads
- `0460901b` refactor: remove unused graph reader paths and payload accessors

Current state: the commands read Newick and FASTA only. The format parsers for MAT (`packages/util-usher-mat`) and PhyloXML (`packages/util-phyloxml`) remain, without graph converters. Tree output goes through `packages/app-output`.

## Decision needed

The background sections of the decision (format landscape and format descriptions) remain valid. The implementation sections and the stated consequence "Users can convert between formats without external tools" do not. Changing a recorded decision needs approval:

- Rewrite the implementation sections to describe the current output path and the parsers without converters
- Or mark the converter design as withdrawn and move the remaining format plans to [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md)

## Related

- [N-io-usher-mat-input-not-wired.md](N-io-usher-mat-input-not-wired.md)
- [N-io-phyloxml-crate-unused.md](N-io-phyloxml-crate-unused.md)
