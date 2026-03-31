# Duplicate edge-collapse implementations in prune and optimize

`collapse_sparse_edge()` in the prune command and `collapse_edge_for_optimize()` in the optimize command implement near-identical edge collapse logic. Both sum branch lengths, compose substitutions via `compose_substitutions()`, and merge indels. The optimize version also handles dense partitions and removes stale partition entries. A shared `topology_cleanup` module would eliminate the duplication and ensure both commands stay consistent.

## Locations

- Prune: [packages/treetime/src/commands/prune/run.rs#L281-L311](../../packages/treetime/src/commands/prune/run.rs#L281-L311) (`collapse_sparse_edge()`) -- sparse partitions only
- Optimize: [packages/treetime/src/commands/optimize/run.rs#L294-L339](../../packages/treetime/src/commands/optimize/run.rs#L294-L339) (`collapse_edge_for_optimize()`) -- sparse + dense, cleans up stale partition entries

Both import `compose_substitutions()` from [packages/treetime/src/seq/mutation.rs#L72-L123](../../packages/treetime/src/seq/mutation.rs#L72-L123). The optimize command already imports `merge_shared_mutation_branches()` from the prune module.

## Related duplication

The deslop audit identified additional structural duplication:

- `write_graph()` duplicated across ancestral, optimize, and prune commands
- GTR initialization blocks (Fitch compression + dummy JC69 + partition creation) duplicated with FIXME comments across optimize and prune

## v0 comparison

v0 has a single `prune_short_branches()` at `treeanc.py:1475:` used from both `optimize_tree` and the timetree pre-loop. Edge collapse is handled by Bio.Phylo's `collapse()` method. v0 does not have shared-mutation merging. The duplication in v1 arose because optimize and prune were implemented independently with different partition type requirements.

## Design decision

Keep `prune` as a separate user-facing command (useful standalone without alignment for prune-by-name). Share the implementation via a common module. The [command-relationships report](../reports/command-relationships/_index.md) describes the ideal hierarchy:

```
topology_cleanup/
  collapse_edge()           -- compose substitutions, sum branch lengths, merge indels
  prune_short_branches()    -- find + collapse zero-length internal edges
  merge_shared_mutations()  -- find + merge sibling pairs sharing substitutions
```

Both `prune` and `optimize` delegate to this module. The `prune` command is a thin CLI wrapper; the optimize loop calls the same functions inside each iteration.

## Scientific background

Edge collapse is a graph contraction operation where an internal node is removed and its children become children of its parent. When the collapsed edge carries substitutions, these must be composed with child edge substitutions using the Markov semigroup property: if the parent edge maps state $A \to G$ and the child edge maps $G \to T$, the net effect is $A \to T$, not the union $\{A \to G, G \to T\}$. Reversions ($A \to G$ then $G \to A$) cancel to no net change. `compose_substitutions()` implements this correctly.

## Related

- [Command-relationships report](../reports/command-relationships/_index.md) -- ideal architecture with shared `topology_cleanup` module
- [Iterative tree refinement: status and recommendations](../reports/iterative-tree-refinement/10-status-and-recommendations.md) -- recommendation R3

## Implementation pointers

The optimize version (`collapse_edge_for_optimize`) is the superset -- it handles both sparse and dense partitions and cleans up stale node/edge entries. Use it as the canonical implementation. The prune version (`collapse_sparse_edge`) is a strict subset.

Module placement: `packages/treetime/src/representation/topology/` or `packages/treetime/src/commands/shared/` are natural candidates. The functions operate on `GraphAncestral` and partition types, so `representation/` aligns better with the layer structure (commands depend on representation, not the reverse).

The refactoring scope includes three related duplications found by the deslop audit: `write_graph()` across commands, and GTR initialization blocks. These can be addressed in the same module extraction pass or as follow-up commits.

Use `ast-grep` for bulk updates of import paths after moving functions.
