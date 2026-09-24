# Timetree aborts when the coalescent counts lineages over an undated node

`treetime timetree` fails on `data/lassa/L/20` after the forward pass leaves some nodes without an inferred time:

```
Error:
   0: Failed to compute coalescent lineage counts
   1: Coalescent lineage count requires an inferred time for every node, but node
      (key=GraphNodeKey(29)) has none. The coalescent model was likely built before node
      times were recomputed for the current tree topology.

Location: packages/treetime/src/coalescent/events.rs:79
```

## Symptom and reproduction

```bash
treetime timetree --tree=data/lassa/L/20/tree.nwk --dates=data/lassa/L/20/metadata.tsv \
  --aln=data/lassa/L/20/aln.fasta.xz --output-all=<dir>
```

The run reaches the abort after warning about empty time distributions:

```
Timetree forward pass: node 'NODE_0000017' has an empty time distribution; no date was assigned.
The messages meeting at this node leave no time with any probability: the dates below it and the
times the rest of the tree implies have disjoint support.
```

This is the `timetree/lassa/L/20/basic` case in `dev/compare-baseline`.

## Impact and scope

- The command produces no output for this dataset, so a user with dates that conflict with the divergence signal gets an abort instead of a result plus a diagnosis.
- The case is excluded from every before/after output comparison, leaving this path without regression coverage.

## Mechanism

Two separate contracts meet here. The forward pass tolerates an empty posterior: when the message from the rest of the tree and the node's own date constraint have disjoint support, it warns and leaves the node undated. `fn collect_tree_events()` [packages/treetime/src/coalescent/events.rs#L78](../../packages/treetime/src/coalescent/events.rs#L78) then requires a time on every node to order the coalescent events, and errors on the first undated one.

The error text attributes the missing time to a stale coalescent model built before a topology change, which is one way to reach this state but not the one this dataset takes: here the time is missing because the posterior was empty, and the message is therefore misleading.

Deciding the intended behavior is part of the fix: either an empty posterior must be resolved to a time (falling back to the constraint or to the parent's implied time), or the coalescent must define what an undated node contributes to the lineage count, or the disjoint-support condition must fail earlier with a diagnosis naming the conflicting dates.

## Related issues

- [H-timetree-mass-sizing-node-times-break-downstream-invariants.md](H-timetree-mass-sizing-node-times-break-downstream-invariants.md)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
