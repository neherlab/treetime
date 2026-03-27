# Prune command applies merge before prune, should be reversed

The prune command runs `merge_shared_mutation_branches()` before `prune_nodes()`. The correct order is prune first, then merge: collapsing short branches creates larger polytomies, which exposes more shared mutations for merging.

## Current behavior

`run_prune()` at [packages/treetime/src/commands/prune/run.rs#L101-L107](../../packages/treetime/src/commands/prune/run.rs#L101-L107):

```rust
if *merge_shared_mutations {
    merge_shared_mutation_branches(&mut graph, &partitions)?;  // merge FIRST
    graph.build()?;
}
prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;  // prune SECOND
```

## Expected behavior

Prune first, then merge:

```rust
prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;
graph.build()?;

if *merge_shared_mutations {
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;
}
```

## Scientific rationale

Consider a tree where an internal node I has two children (binary resolution), and the edge to I has branch length below the `--prune-short` threshold:

```
Before prune:
  P
  +-- I (bl=1e-8)
  |   +-- A (subs: G100T)
  |   +-- B (subs: G100T)
  +-- C (subs: G100T)
  +-- D (subs: T200A)

After prune (I collapsed):
  P  (polytomy of degree 4)
  +-- A (subs: G100T)
  +-- B (subs: G100T)
  +-- C (subs: G100T)
  +-- D (subs: T200A)

After merge (A, B, C grouped):
  P
  +-- N (subs: G100T)
  |   +-- A
  |   +-- B
  |   +-- C
  +-- D (subs: T200A)
```

Without pruning first, the node I separates A and B from C. The merge step only sees A and B as siblings sharing G100T, and misses that C also shares it. Pruning I first creates a 4-way polytomy where the merge correctly groups A, B, and C together.

Tree builders produce arbitrary binary resolutions of polytomies. Short branches between such arbitrary resolutions are the primary source of hidden polytomies. Removing them first maximizes the number of siblings available for shared-mutation merging.

## Fix

Swap the order of `merge_shared_mutation_branches` and `prune_nodes` in `run_prune()`. The `graph.build()` call after pruning is needed to update the graph's internal index before the merge step iterates over polytomy nodes.

When both operations are integrated into the optimize loop (see [M-optimize-no-topology-cleanup-in-loop.md](M-optimize-no-topology-cleanup-in-loop.md)), the same prune-then-merge order must be maintained in each iteration.

## Dependencies

- None. This is a standalone fix to the prune command.
- If `--prune-short` edges carry mutations, the composition bug in [M-prune-collapse-uses-union-not-composition.md](M-prune-collapse-uses-union-not-composition.md) affects the mutation lists passed to the merge step. Fixing composition first ensures the merge step sees correct substitution data.

## Impact

Low for `--prune-empty` (zero-mutation edges produce no polytomy enlargement benefit). Medium for `--prune-short` with alignments where tree builders introduce many short arbitrary resolutions (flu, SARS-CoV-2 datasets with many identical or near-identical sequences).
