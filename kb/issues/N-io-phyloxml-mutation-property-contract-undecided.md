# PhyloXML mutation property contract is undecided

## Problem

PhyloXML does not define a native mutation element. Its generic `Property` type carries typed free text, a `prefix:value` reference key, an `applies_to` target, and an optional `id_ref` that attaches the property to a specific XML element [packages/util-phyloxml/schemas/phyloxml.xsd#L385-L418](../../packages/util-phyloxml/schemas/phyloxml.xsd#L385-L418).

TreeTime currently writes substitutions as `gene:A123T` and indels as `gene:ins|del:start:SEQ` under ad hoc `treetime:mutation` and `treetime:indel` references [packages/treetime-io/src/tree_ir/phyloxml.rs#L217-L260](../../packages/treetime-io/src/tree_ir/phyloxml.rs#L217-L260). These values round-trip through TreeTime, but the reference vocabulary is undocumented and the two event kinds use unrelated grammars. External consumers have no published contract for interpreting them.

The design space was already identified in [kb/proposals/phyloxml-treetime-property-namespace.md](../proposals/phyloxml-treetime-property-namespace.md), but its actionable decisions were not represented by a known issue.

## Potential solutions

### A1. Property reference vocabulary

- **Keep the undeclared pseudo-prefix:** retain TreeTime round trips but leave properties non-self-describing.
- **Publish a stable TreeTime reference vocabulary:** define every `treetime:*` key and associate the prefix with durable documentation. XML namespace declarations do not bind ordinary `Property.ref` attribute values.

### A2. Mutation values

- **Keep separate bespoke substitution and indel strings:** preserves current round trips but duplicates the application mutation grammar.
- **Use canonical atomic mutation strings:** emit one property per substitution, inserted site, or deleted site; include the gene track in a defined wrapper.
- **Use one self-contained structured value per event:** encode event kind, track, position/range, and sequence in a published string or structured-data grammar. `id_ref` cannot group sibling properties because it identifies the XML element to which a property applies.

### A3. Round-trip ownership

- **TreeTime-private extension:** optimize only for TreeTime read/write symmetry.
- **Published interoperability contract:** document reference keys, value grammar, coordinate base, repeatability, and unknown-property behavior for other consumers.

## Recommendation

Publish a stable TreeTime reference vocabulary and interoperability contract. Reuse canonical atomic mutation strings for substitutions and aligned per-site gap changes, with the track wrapper defined by that contract. If grouped indel event identity is required, use one self-contained value per event; do not misuse `id_ref` as a grouping identifier.

## Ticket readiness

No implementation ticket is ready. Reference-vocabulary identity, atomic versus grouped indel representation, and external compatibility behavior require explicit approval. This issue remains the tracking item until those axes are decided.

## Related items

- [kb/proposals/phyloxml-treetime-property-namespace.md](../proposals/phyloxml-treetime-property-namespace.md)
- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
