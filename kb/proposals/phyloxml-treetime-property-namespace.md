# A dedicated PhyloXML property namespace for TreeTime data

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

These keys round-trip within TreeTime but are not a registered or documented namespace, so other PhyloXML consumers cannot interpret them, and the encoding of compound values (mutation strings, indel tuples) is bespoke.

## Open questions

- **Namespace registration**: PhyloXML `<property>` `ref` values are conventionally namespaced as `NS:key` where `NS` is declared in the document. Should TreeTime declare an `xmlns:treetime` namespace URI and register the property set, so the document is self-describing?
- **Compound value encoding**: mutations and indels are currently packed into delimited strings. PhyloXML has no list type, but a structured alternative (one property per field with shared `id_ref`, or repurposing `<sequence>`/`<binary_characters>`) could be more interoperable.
- **Alignment with augur/auspice**: auspice represents mutations as `gene -> [A123T, ...]`. Should the PhyloXML encoding mirror that grouping to ease cross-format conversion?

## Non-goals

This proposal does not change the round-trip behavior within TreeTime, which is already correct. It concerns interoperability and self-description for external PhyloXML consumers.
