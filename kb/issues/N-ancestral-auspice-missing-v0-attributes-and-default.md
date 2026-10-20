# Ancestral Auspice output lacks v0 node attributes and is not a default output

> [!IMPORTANT]
> **Decision required.** The Auspice payload (A2) is unapproved. Matching v0 adds `confidence` and `bad_branch`; a subset without them is an intentional parity divergence. See "Decision axes".

V0's `ancestral` command always writes `auspice_tree.json` through `def export_sequences_and_tree()`. V1 writes Auspice output for every reconstruction method (Fitch parsimony, sparse marginal, dense marginal) directly from the graph via `fn auspice_tree()` [packages/app-output/src/auspice.rs#L29](../../packages/app-output/src/auspice.rs#L29). `fn write_ancestral_trees()` [packages/app-commands/src/commands/shared/ancestral_trees.rs#L20-L69](../../packages/app-commands/src/commands/shared/ancestral_trees.rs#L20-L69) passes the reconstructed sequences and branch mutations of each partition to the tree writers.

V0's ancestral Auspice JSON contains `node_attrs.div` (cumulative `mutation_length`), `branch_attrs.mutations.nuc` (per-branch mutations), `node_attrs.confidence` (pseudo-bootstrap from the number of mutations), `meta.genome_annotations.nuc`, and `node_attrs.bad_branch`. It contains no dates because `timetree=false`.

V1 writes `div`, the nucleotide branch mutations, and `meta.genome_annotations.nuc`. Two parts remain:

- **Missing node attributes**: v1 writes no `confidence` attribute, and writes `bad_branch` only when the graph carries dates [packages/app-output/src/auspice.rs#L227-L230](../../packages/app-output/src/auspice.rs#L227-L230), which `ancestral` never does
- **Not in the defaults**: `fn CommandKind.default_outputs()` [packages/app-output/src/output_plan.rs#L597-L616](../../packages/app-output/src/output_plan.rs#L597-L616) adds Auspice to the tree defaults only for `timetree`, so `--output-all` without a selection writes no Auspice file for `ancestral`. The file is written only when the user selects it

This is a standalone visualization convenience for viewing ancestral reconstruction results in auspice without going through `augur export v2`. The augur pipeline path is covered by node data JSON output.

## Decision axes

### A2. Auspice payload

- O1. Match the v0 ancestral payload: cumulative `div`, nucleotide branch mutations, confidence, nucleotide genome annotation, and `bad_branch`.
- O2. Emit the current subset and document the omitted inference fields as an intentional divergence.

**Recommendation:** O1, the project default of reference parity. Shared metadata construction should be reused with timetree where the field semantics are identical.

O1 adds the `bad_branch` filter, which then becomes the default coloring in place of the most recurrent site (`fn genotype_color_by()` [packages/app-output/src/auspice.rs#L127](../../packages/app-output/src/auspice.rs#L127)). Every ancestral branch has the value `No`, so the tree opens in one color.

## Evidence

- `def create_auspice_json()` [packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341)
- `def export_sequences_and_tree()` Auspice dispatch [packages/legacy/treetime/treetime/CLI_io.py#L220](../../packages/legacy/treetime/treetime/CLI_io.py#L220)
- `def ancestral_reconstruction()` output dispatch [packages/legacy/treetime/treetime/wrappers.py#L641-L648](../../packages/legacy/treetime/treetime/wrappers.py#L641-L648)
- `fn write_ancestral_trees()` [packages/app-commands/src/commands/shared/ancestral_trees.rs#L20-L69](../../packages/app-commands/src/commands/shared/ancestral_trees.rs#L20-L69)
- `fn CommandKind.available_tree_outputs()` [packages/app-output/src/output_plan.rs#L618-L628](../../packages/app-output/src/output_plan.rs#L618-L628)
- `fn CommandKind.default_outputs()` [packages/app-output/src/output_plan.rs#L597-L616](../../packages/app-output/src/output_plan.rs#L597-L616)

## Validation

- Exercise Fitch, sparse marginal, and dense marginal ancestral runs with `--output-all`, and assert that the Auspice file is written.
- Parse and schema-validate the Auspice artifact, then compare semantic fields with the v0 oracle.

## Related

- [N-io-auspice-output-missing-augur-golden-master.md](N-io-auspice-output-missing-augur-golden-master.md) - golden master of the whole Auspice document
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md) - node data JSON shares the same data sources (partition mutations, annotations)
