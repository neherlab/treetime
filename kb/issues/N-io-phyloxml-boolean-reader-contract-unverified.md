# PhyloXML boolean reader contract is unverified after the tree I/O refactor

The removed `tree_ir` PhyloXML reader parsed recognized boolean properties with `value == "true"`, which coerced the XML Schema value `1` and every invalid lexical form to false. The current generic reader in [`packages/treetime-io/src/phyloxml.rs`](../../packages/treetime-io/src/phyloxml.rs) delegates node conversion through `PhyloxmlToGraph`; current source inspection does not establish a command-level conversion for `treetime:bad_branch` or `treetime:date_inferred`.

The historical coercion therefore does not establish current incorrect behavior. The remaining problem is an unverified round-trip and validation contract for recognized properties.

## Evidence required

- Identify every current `PhyloxmlToGraph` implementation and whether the recognized TreeTime properties are reachable through a public command.
- Verify the supported XML Schema boolean lexical forms: `true`, `false`, `1`, and `0`.
- Determine whether invalid recognized values must be rejected and which conversion boundary can preserve property `ref` context in errors.
- Verify `bad_branch` and `date_inferred` round-trip semantics before prescribing a parser helper.

## Ticket readiness

No implementation ticket is ready. The removed `tree_ir` paths do not identify the current conversion boundary.
