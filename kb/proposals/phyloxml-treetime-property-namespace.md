# A dedicated PhyloXML property vocabulary for TreeTime data

## Motivation

PhyloXML 1.20 models a fixed set of biological elements (name, branch length, taxonomy, sequence, date, confidence, distribution). TreeTime carries data that has no native PhyloXML element: cumulative divergence, discrete trait assignments, per-branch nucleotide and amino-acid substitutions, indels, the relaxed-clock branch-rate multiplier (gamma), and the bad-branch (temporal outlier) flag.

The PhyloXML schema provides a generic `<property>` element for exactly this purpose. Each property carries a free-form `ref` key, a `datatype` (XSD type), an `applies_to` qualifier (`node`, `parent_branch`, `clade`, `phylogeny`, ...), and a string value. The TreeIR PhyloXML adapter projects TreeTime data into these properties.

## Current state

The adapter (`packages/treetime-io/src/tree_ir/phyloxml.rs`) uses ad-hoc `ref` keys under a `treetime:` pseudo-prefix:

- `treetime:divergence` (xsd:double, node)
- `treetime:bad_branch` (xsd:boolean, node)
- `treetime:date_inferred` (xsd:boolean, node)
- `treetime:gamma` (xsd:double, parent_branch)
- `treetime:mutation` (xsd:string, parent_branch), value `gene:A123T`
- `treetime:indel` (xsd:string, parent_branch), value `gene:ins|del:start:SEQ`
- `treetime:trait:<attr>` (xsd:string, node)

Scalar keys have reader and writer paths within TreeTime, but the mutation and indel vocabulary, event grouping, and round-trip ownership contract remain undecided. The keys are not a registered or documented reference vocabulary, so other PhyloXML consumers cannot interpret them, and the encoding of compound values is bespoke. `Property.ref` is a constrained attribute value rather than an XML QName, so declaring `xmlns:treetime` would not define its semantics [packages/util-phyloxml/schemas/phyloxml.xsd#L385-L408](../../packages/util-phyloxml/schemas/phyloxml.xsd#L385-L408).

## Open questions

- **Reference vocabulary**: PhyloXML `<property>` `ref` values use `prefix:value` syntax. What durable documentation URI and versioning contract should define the `treetime:*` vocabulary?
- **Compound value encoding**: mutations and indels are currently packed into delimited strings. PhyloXML has no list type. Canonical atomic values, one self-contained structured value per event, or a schema extension could be more interoperable. `id_ref` targets an XML element and cannot group sibling properties into an event.
- **Alignment with augur/auspice**: auspice represents mutations as `gene -> [A123T, ...]`. Should the PhyloXML encoding mirror that grouping to ease cross-format conversion?

## Non-goals

This proposal does not change native PhyloXML elements or scalar TreeTime properties. It defines the unresolved interoperability, self-description, and round-trip contract for typed mutations and indels.

## Related issues

- [kb/issues/N-io-phyloxml-mutation-property-contract-undecided.md](../issues/N-io-phyloxml-mutation-property-contract-undecided.md)
