# Ancestral Auspice output lacks v0 node attributes and is not a default output

> [!IMPORTANT]
> **Decision required.** The Auspice payload (A2) is unapproved. Matching v0 adds `confidence` and `bad_branch`; a subset without them is an intentional parity divergence. See "Decision axes".

V0's `ancestral` command always writes `auspice_tree.json` through `def export_sequences_and_tree()`. V1 writes Auspice output for every reconstruction method (Fitch parsimony, sparse marginal, dense marginal) directly from the graph via `fn auspice_tree()` [packages/app-output/src/auspice.rs#L32](../../packages/app-output/src/auspice.rs#L32). `fn write_ancestral_trees()` [packages/app-commands/src/commands/shared/ancestral_trees.rs#L21-L71](../../packages/app-commands/src/commands/shared/ancestral_trees.rs#L21-L71) passes the reconstructed sequences and branch mutations of each partition to the tree writers.

V0's ancestral Auspice JSON contains `node_attrs.div` (cumulative `mutation_length`), `branch_attrs.mutations.nuc` (per-branch mutations), `node_attrs.confidence` (pseudo-bootstrap from the number of mutations), `meta.genome_annotations.nuc`, and `node_attrs.bad_branch`. It contains no dates because `timetree=false`.

V1 writes `div`, the nucleotide branch mutations, and `meta.genome_annotations.nuc`. Two parts remain:

- **Missing node attributes**: v1 writes no `confidence` attribute, and writes `bad_branch` only when the graph carries dates [packages/app-output/src/auspice.rs#L254-L257](../../packages/app-output/src/auspice.rs#L254-L257), which `ancestral` never does
- **Not in the defaults**: `fn CommandKind.outputs()` [packages/app-output/src/command_outputs.rs#L42-L56](../../packages/app-output/src/command_outputs.rs#L42-L56) marks Auspice `DefaultRule::Default` only for `timetree` (`OnRequest` for `ancestral`), so `--output-all` without a selection writes no Auspice file for `ancestral`. The file is written only when the user selects it

This is a standalone visualization convenience for viewing ancestral reconstruction results in auspice without going through `augur export v2`. The augur pipeline path is covered by node data JSON output.

## Decision axes

### A2. Auspice payload

- O1. Match the v0 ancestral payload: cumulative `div`, nucleotide branch mutations, confidence, nucleotide genome annotation, and `bad_branch`.
- O2. Emit the current subset and document the omitted inference fields as an intentional divergence.

**Recommendation:** O1, the project default of reference parity. Shared metadata construction should be reused with timetree where the field semantics are identical.

O1 adds the `bad_branch` filter, which then becomes the default coloring in place of the most recurrent site (`fn genotype_color_by()` [packages/app-output/src/auspice.rs#L148](../../packages/app-output/src/auspice.rs#L148)). Every ancestral branch has the value `No`, so the tree opens in one color.

## Evidence

- `def create_auspice_json()` [packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341)
- `def export_sequences_and_tree()` Auspice dispatch [packages/legacy/treetime/treetime/CLI_io.py#L220](../../packages/legacy/treetime/treetime/CLI_io.py#L220)
- `def ancestral_reconstruction()` output dispatch [packages/legacy/treetime/treetime/wrappers.py#L641-L648](../../packages/legacy/treetime/treetime/wrappers.py#L641-L648)
- `fn write_ancestral_trees()` [packages/app-commands/src/commands/shared/ancestral_trees.rs#L21-L71](../../packages/app-commands/src/commands/shared/ancestral_trees.rs#L21-L71)
- `fn CommandKind.outputs()` [packages/app-output/src/command_outputs.rs#L42-L56](../../packages/app-output/src/command_outputs.rs#L42-L56)
- `fn CommandKind.default_outputs()` [packages/app-output/src/command_outputs.rs#L26](../../packages/app-output/src/command_outputs.rs#L26)

## Validation

- Exercise Fitch, sparse marginal, and dense marginal ancestral runs with `--output-all`, and assert that the Auspice file is written.
- Parse and schema-validate the Auspice artifact, then compare semantic fields with the v0 oracle.

## Related

- [N-io-auspice-output-missing-augur-golden-master.md](N-io-auspice-output-missing-augur-golden-master.md) - golden master of the whole Auspice document
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md) - node data JSON shares the same data sources (partition mutations, annotations)

## Proposed implementation plan

The payload (A2) stays open; the default-output step needs no payload answer.

1. **Default output** (both options): change the `ancestral` entry of `OutputSelection::Auspice` in `fn CommandKind.outputs()` [packages/app-output/src/command_outputs.rs#L47](../../packages/app-output/src/command_outputs.rs#L47) from `DefaultRule::OnRequest` to `DefaultRule::Default`, as v0 always writes `auspice_tree.json`. Update the `--help` text and docs that list the default outputs, then `just gen`
2. **O1**:
   - `confidence`: compute the v0 pseudo-bootstrap value from the number of mutations per branch, with the formula of `def create_auspice_json()` [packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341), as a new node attribute in `fn auspice_tree()` [packages/app-output/src/auspice.rs#L32](../../packages/app-output/src/auspice.rs#L32). The field must not collide with the date `confidence` of `num_date` (line 252) or the trait confidence (line 344)
   - `bad_branch`: write `No` for every node when the graph has no dates ([packages/app-output/src/auspice.rs#L254-L257](../../packages/app-output/src/auspice.rs#L254-L257)), and check the coloring order of `fn genotype_color_by()` (line 148) against v0
3. **O2**: write a `kb/decisions/` entry after approval and keep the code

Tests (oracle: v0 `treetime ancestral` output on the same input):

- Unit test in `packages/app-output/src/__tests__/` for the confidence value on a small tree with known mutation counts, expected values computed by hand from the v0 formula
- Default-output test: `CommandKind::Ancestral.default_outputs()` contains `OutputSelection::Auspice`
- Validation run: Fitch, sparse marginal, and dense marginal `just r treetime ancestral --output-all` on `flu/h3n2/20`; compare `div`, mutations, `confidence`, and `bad_branch` with the v0 file. Then `just check`, `just test-rs`, `just generated-check`
