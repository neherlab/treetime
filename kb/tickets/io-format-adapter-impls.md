# Wire per-command TreeIR projections for PhyloXML, Auspice, and UShER MAT

## Status

The format-neutral TreeIR foundation is implemented and tested. Remaining work is the per-command projection of analysis results into the IR graph and enabling the formats in the availability matrix.

Done:

- `packages/treetime-io/src/tree_ir/`: `TreeIrNode`/`TreeIrEdge`/`TreeIrData`, mutation/indel model (`TreeIrSub`/`TreeIrIndel` over `AsciiChar`), and bidirectional PhyloXML, Auspice v2, and UShER MAT adapters, with round-trip and field-mapping tests.
- `write_tree_outputs` takes an optional `TreeIrGraph`; PhyloXML/Auspice/MAT route through it, Newick/Nexus/Dot/graph-JSON stay on the domain graph. The old `AuspiceWriter` trait and stubs are removed.
- Timetree writes Auspice through the IR (`commands/timetree/output/ir.rs`), replacing `build_timetree_auspice`/`TimetreeAuspiceWriter`. Verified by unit tests and an end-to-end smoke run.

## Remaining task

For ancestral, optimize, mugration, clock, and prune: build a `TreeIrGraph` from the command result + partition data at write time, pass it to `write_tree_outputs`, and enable the command's formats in `CommandKind::available_tree_outputs()` per the matrix below. Extend timetree's projection with per-branch mutations to enable its PhyloXML/MAT outputs.

| format           | ancestral | optimize | timetree | mugration | clock | prune |
| ---------------- | --------- | -------- | -------- | --------- | ----- | ----- |
| phyloxml (+json) | yes       | yes      | yes      | yes       | yes   | yes   |
| auspice          | yes       | no       | yes      | yes       | no    | no    |
| mat-pb (+json)   | yes       | yes      | yes      | no        | no    | no    |

Use `commands/timetree/output/ir.rs` as the projection pattern.

## Projection data sources

- **ancestral** (`commands/ancestral/run.rs`): name, desc, confidence. Nucleotide mutations via `MutationCommentProvider` / `ml_subs()`/`fitch_subs()` (`partition/traits.rs`); amino-acid mutations via `AaNodeData` (`commands/ancestral/aa_node_data.rs`); indels via `SparseEdgePartition.indels`; root sequence from partition; divergence from `compute_edge_mutation_counts` (`seq/div.rs`).
- **optimize** (`commands/optimize/run.rs`): name, branch lengths, mutations (same extraction as ancestral).
- **mugration** (`commands/mugration/run.rs`): discrete trait values via `DiscreteCommentProvider` (`partition/marginal_discrete.rs`) -> `TreeIrTrait` (value, confidence map, entropy). No mutations.
- **clock** (`commands/clock/run.rs`): name, `node.time`, `node.div`, bad-branch/outlier. PhyloXML only.
- **prune** (`commands/prune/run.rs`): name, branch lengths. PhyloXML only.

When mapping mutations, convert `Sub` -> `TreeIrSub` (`AsciiChar`); `treetime-io` must not depend on `treetime`.

## Tests to add

- Round-trip and field-mapping per command (extend the `treetime-io` IR tests where format-level).
- Snapshot tests: small fixed trees -> PhyloXML and Auspice golden strings.
- Golden master: ancestral/mugration Auspice and MAT against augur/UShER where an oracle exists.
- Edge cases: indels to UShER (warning + drop), 100+ mutations/branch, empty/single-node tree, no-name/no-branch-length nodes, large-tree smoke.
- Smoke tests: ancestral/timetree with `--output-tree-phyloxml`, `--output-tree-auspice`, `--output-tree-mat-pb`.

## Related issues

- Source: [kb/issues/N-io-write-graph-files-missing-formats.md](../issues/N-io-write-graph-files-missing-formats.md)
