# Embedding substitution and clock models in tree output formats

## Motivation

TreeTime infers a GTR substitution model and a molecular-clock model alongside the tree. The analysis commands write these to sidecar files (`--output-gtr`, `--output-clock-model`). Some tree formats can carry model parameters inline, which would make a single output file self-contained for downstream tools.

The graph-backed output path (PhyloXML, Auspice v2, UShER MAT) currently carries only per-node and per-branch data plus minimal graph-level metadata. It does not embed the GTR or clock model.

## Options

- **PhyloXML**: no native model element; would require `<property>` elements at `applies_to=phylogeny` scope (e.g. a serialized GTR matrix and clock rate). Verbose and non-standard.
- **Auspice v2**: the schema has no substitution-model field. Clock information is implicitly present through `num_date`. augur keeps the model out of the auspice JSON; embedding it would diverge from augur.
- **UShER MAT**: the protobuf has no model field. Out of scope for the format.

## Recommendation (open)

Keep models in sidecar files (current behavior) until a concrete downstream need justifies embedding. If embedding is pursued, prefer PhyloXML phylogeny-scope properties and a documented schema, and never diverge from augur's auspice output for the auspice path.

## Non-goals

This proposal does not change current sidecar output. It records that inline model embedding is absent from graph-backed tree output.
