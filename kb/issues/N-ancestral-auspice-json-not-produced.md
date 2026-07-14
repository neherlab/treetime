# Ancestral Auspice output is incomplete and method-dependent

V0's `ancestral` command produces `auspice_tree.json` through `def export_sequences_and_tree()`. V1 marginal reconstruction constructs TreeIR and can produce Auspice output. Parsimony advertises TreeIR-backed formats without constructing the required projection, so a valid selection can fail after earlier outputs have been written.

V0's ancestral Auspice JSON contains `node_attrs.div` (cumulative `mutation_length`), `branch_attrs.mutations.nuc` (per-branch mutations), `node_attrs.confidence` (pseudo-bootstrap), `meta.genome_annotations.nuc`, and `node_attrs.bad_branch`. It contains no dates because `timetree=false`.

This is a standalone visualization convenience for viewing ancestral reconstruction results in auspice without going through `augur export v2`. The augur pipeline path is covered by node data JSON output.

## Decision axes

### A1. Method availability

- O1. Construct the same typed TreeIR projection for Fitch and marginal partitions. This matches v0's method-independent output capability and keeps the advertised format matrix accurate.
- O2. Restrict Auspice and other TreeIR-backed formats to marginal reconstruction and reject a parsimony selection during output planning, before any file is written. This is internally consistent but is a parity divergence requiring explicit approval.

**Recommendation:** O1. Fitch partitions already contain resolved node states and edge mutations, so the writer should receive those values through the same typed projection contract.

### A2. Auspice payload

- O1. Match the v0 ancestral payload: cumulative `div`, nucleotide branch mutations, confidence, nucleotide genome annotation, and `bad_branch`.
- O2. Emit a schema-valid subset and document the omitted inference fields as an intentional divergence.

**Recommendation:** O1, the project default of reference parity. Shared metadata construction should be reused with timetree where the field semantics are identical.

### A3. Failure atomicity

- O1. Validate format availability and construct every requested projection before opening output files.
- O2. Keep sequential writer dispatch and remove files after a later projection failure. Cleanup can itself fail and exposes partial output during execution.

**Recommendation:** O1. Treat output planning and projection as a fallible preparation step, then perform writes only after every requested artifact is known to be constructible.

## Evidence

- `def create_auspice_json()` [packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341)
- `def export_sequences_and_tree()` Auspice dispatch [packages/legacy/treetime/treetime/CLI_io.py#L220](../../packages/legacy/treetime/treetime/CLI_io.py#L220)
- `def ancestral_reconstruction()` output dispatch [packages/legacy/treetime/treetime/wrappers.py#L641-L648](../../packages/legacy/treetime/treetime/wrappers.py#L641-L648)
- `fn write_tree_for_partition()` [packages/treetime/src/commands/ancestral/run.rs#L235-L270](../../packages/treetime/src/commands/ancestral/run.rs#L235-L270)
- `pub fn write_tree_outputs()` and `fn require_ir()` [packages/treetime-io/src/graph.rs#L61-L126](../../packages/treetime-io/src/graph.rs#L61-L126)
- `fn CommandKind.available_tree_outputs()` [packages/treetime/src/commands/shared/output.rs#L297-L322](../../packages/treetime/src/commands/shared/output.rs#L297-L322)

## Validation

- Exercise Fitch, sparse marginal, and dense marginal ancestral runs for every advertised format.
- Parse and schema-validate the Auspice artifact, then compare semantic fields with the v0 oracle.
- Inject a projection failure and assert that no requested output path has been created or modified.

## Ticket readiness

No implementation ticket is ready. O1 is the reference-parity recommendation, but implementing it through TreeIR depends on approval of the shared TreeIR boundary. O2 is an intentional parity divergence and also requires approval. The independent topology-ordering defect has its own ready ticket.

## Related

- [M-timetree-tree-output-inference-metadata-incomplete.md](M-timetree-tree-output-inference-metadata-incomplete.md) - shared mutation, confidence, and annotation contract
- [M-io-tree-backed-output-order-inconsistent.md](M-io-tree-backed-output-order-inconsistent.md) - shared topology ordering
- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md) - typed mutation projection
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md) - node data JSON shares the same data sources (partition mutations, annotations)
