# Move `merge_shared_mutation_branches()` to shared topology_cleanup module

`merge_shared_mutation_branches()` is defined in `commands/prune/run.rs` and imported by `commands/optimize/run.rs` via a cross-command `pub(crate) use` path. The edge-collapse counterpart has already been extracted to `representation/algo/topology_cleanup/collapse`; `merge_shared_mutation_branches` should follow for architectural symmetry and to remove the cross-command import.

## Locations

- Definition: [packages/treetime/src/commands/prune/run.rs#L293-L313](../../packages/treetime/src/commands/prune/run.rs#L293-L313) (plus the helper functions `find_polytomy_nodes`, `merge_single_polytomy`, `collect_child_edge_keys`, `find_best_shared_mutation_pair`, `compute_shared_subs_across_partitions`, and `merge_sibling_pair`)
- Optimize-side usage: [packages/treetime/src/commands/optimize/run.rs](../../packages/treetime/src/commands/optimize/run.rs) -- `use crate::commands::prune::run::merge_shared_mutation_branches;`

## Proposed placement

`packages/treetime/src/representation/algo/topology_cleanup/merge_shared_mutations.rs` -- alongside `collapse.rs`. The `prune` command reduces to a thin CLI wrapper that imports the shared function; the optimize loop imports it directly.

## Related

- [Command-relationships report](../reports/command-relationships/_index.md) -- describes the ideal `topology_cleanup` module containing `collapse_edge()`, `prune_short_branches()`, and `merge_shared_mutations()`.
- [Iterative tree refinement R3](../reports/iterative-tree-refinement/10-status-and-recommendations.md) -- shared-module extraction recommendation (collapse done, merge pending).
- `packages/treetime/src/representation/algo/topology_cleanup/collapse.rs` -- shared collapse_edge implementation to follow.
