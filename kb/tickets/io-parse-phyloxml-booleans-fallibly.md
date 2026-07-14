# Parse PhyloXML booleans fallibly

Implement the XML Schema boolean lexical contract for recognized PhyloXML properties.

## Required changes

- Add one private fallible parser that maps `true` and `1` to true and maps `false` and `0` to false.
- Reject every other lexical value with the property `ref` and original value in the error.
- Use the parser for `REF_BAD_BRANCH` and `REF_DATE_INFERRED` in `clade_to_graph_components()` [packages/treetime-io/src/tree_ir/phyloxml.rs#L165](../../packages/treetime-io/src/tree_ir/phyloxml.rs#L165).

## Validation

- Parameterize all four legal forms across both recognized boolean properties.
- Cover invalid case variants, surrounding whitespace, an empty value, and unrelated text.
- Assert complete parsed nodes for valid documents and contextual error content for invalid documents.
- Confirm valid documents round-trip without changing either boolean.
- Run the full lint and test suite.

## Related issues

- Source: [kb/issues/M-io-phyloxml-booleans-silently-coerced.md](../issues/M-io-phyloxml-booleans-silently-coerced.md)
