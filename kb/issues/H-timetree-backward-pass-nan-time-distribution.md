# Timetree aborts on lassa/L/20 with a NaN time distribution

`treetime timetree` fails on `data/lassa/L/20` in the first backward pass, when it normalizes the combined time distribution of an internal node:

```
Error:
   0: When normalizing the time distribution of node 37
   1: Cannot normalize a distribution on [314.4969722119549, 1976.1488159980956]: its negative
      log-likelihood values contain NaN

Location: packages/treetime-distribution/src/distribution_core/distribution.rs:215
```

## Symptom and reproduction

```bash
treetime timetree --tree=data/lassa/L/20/tree.nwk --dates=data/lassa/L/20/metadata.tsv \
  --aln=data/lassa/L/20/aln.fasta.xz --output-all=<dir>
```

This is the `timetree/lassa/L/20/basic` case of `dev/smoke`.

Whether the dates of this dataset also have disjoint support is unknown until the `NaN` is removed.

## Impact and scope

- The command produces no output for this dataset.
- The case is excluded from every before/after output comparison, leaving this path without regression coverage.

## Mechanism

The combined backward message at the node is the product of the child messages, the coalescent prior, and the date constraint, so the `NaN` comes from one of those factors or from their product. In negative-log space a product is a sum, and `+inf + -inf` is `NaN`; a message with `-inf` or `NaN` values would produce this result. The source has not been traced.

When the `NaN` is removed, two contracts still meet here. The forward pass tolerates an empty posterior: when the message from the rest of the tree and the node's own date constraint have disjoint support, it warns and leaves the node undated. `fn collect_tree_events()` [packages/treetime/src/coalescent/events.rs#L78](../../packages/treetime/src/coalescent/events.rs#L78) then requires a time on every node to order the coalescent events, and errors on the first undated one. Its error text blames a stale coalescent model, which is misleading when the posterior is empty.

Deciding the intended behavior for that case is part of the fix: either an empty posterior must be resolved to a time (falling back to the constraint or to the parent's implied time), or the coalescent must define what an undated node contributes to the lineage count, or the disjoint-support condition must fail earlier with a diagnosis naming the conflicting dates.

## Related issues

- [H-timetree-mass-sizing-node-times-break-downstream-invariants.md](H-timetree-mass-sizing-node-times-break-downstream-invariants.md)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md)
