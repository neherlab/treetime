# Per-command projections missing for PhyloXML, Auspice, and UShER MAT

The format-neutral TreeIR output path is implemented in `treetime-io`: a `Graph<TreeIrNode, TreeIrEdge, TreeIrData>` with bidirectional PhyloXML, Auspice v2, and UShER MAT adapters, and `write_tree_outputs` routes those three formats through an optional IR graph. Timetree builds an IR graph and writes Auspice through this path (replacing the old `TimetreeAuspiceWriter`).

What remains is the per-command projection of analysis results into the IR graph for the other commands, and enabling the corresponding formats in `CommandKind::available_tree_outputs()` per the availability matrix.

## Remaining work

Each command must build a TreeIR graph from its result + partition data at write time and pass it to `write_tree_outputs`, then enable its formats in `available_tree_outputs()`:

| format           | ancestral | optimize | timetree | mugration | clock | prune |
| ---------------- | --------- | -------- | -------- | --------- | ----- | ----- |
| phyloxml (+json) | todo      | todo     | partial  | todo      | todo  | todo  |
| auspice          | todo      | n/a      | done     | todo      | n/a   | n/a   |
| mat-pb (+json)   | todo      | todo     | todo     | n/a       | n/a   | n/a   |

Timetree Auspice is `done` (div/date/bad-branch via IR). Timetree PhyloXML/MAT are `partial`: the IR projection carries div/date/bad-branch but not per-branch mutations, so those formats are not yet enabled for timetree.

Projection data sources to wire (per `kb/tickets/io-format-adapter-impls.md`):

- Nucleotide substitutions per edge: `MutationCommentProvider` / `ml_subs()`/`fitch_subs()` (`partition/traits.rs`), mapped from `Sub` to `TreeIrSub` (`AsciiChar`).
- Amino-acid substitutions: `AaNodeData` (`commands/ancestral/aa_node_data.rs`).
- Indels: `SparseEdgePartition.indels` -> `TreeIrIndel`.
- Divergence: `compute_edge_mutation_counts` (`seq/div.rs`) for `divergence-units=mutations`, branch length otherwise.
- Discrete traits (mugration): `DiscreteCommentProvider` (`partition/marginal_discrete.rs`) -> `TreeIrTrait`.
- Root reference sequence: partition -> `TreeIrData.root_sequence`.

## Locations

- IR types and adapters: `packages/treetime-io/src/tree_ir/`
- Output dispatch: `packages/treetime-io/src/graph.rs` (`write_tree_outputs`, `TreeIrGraph`)
- Per-command wiring: `packages/treetime/src/commands/{ancestral,optimize,timetree,mugration,clock,prune}/run.rs`
- Format availability: `packages/treetime/src/commands/shared/output.rs` (`available_tree_outputs`)
- Timetree IR projection (reference pattern): `packages/treetime/src/commands/timetree/output/ir.rs`

## Related

- [kb/tickets/io-format-adapter-impls.md](../tickets/io-format-adapter-impls.md) -- implementation instructions
- [kb/proposals/phyloxml-treetime-property-namespace.md](../proposals/phyloxml-treetime-property-namespace.md)
- [kb/proposals/tree-format-model-embedding.md](../proposals/tree-format-model-embedding.md)
