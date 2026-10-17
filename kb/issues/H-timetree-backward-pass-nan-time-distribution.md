# Timetree aborts with a NaN time distribution

> [!WARNING]
> **Needs review.** Commit `193628ba8d` re-peaks the node's combined distribution instead of normalizing it. The same `NaN` is now reported as "When re-peaking the time distribution of node 37: Cannot re-peak a distribution on [...]: its values contain NaN", from `fn DistributionFunction.repeak()` in [packages/treetime-distribution/src/distribution_core/function.rs#L340](../../packages/treetime-distribution/src/distribution_core/function.rs#L340). The quoted error and location below predate it.

> [!IMPORTANT]
> **Investigation required.** The source of the `NaN` is not traced. On one reproduction below, find the first factor of the combined message (child message, coalescent prior, date constraint) that holds `NaN` or `-inf`, and the operation that produces it.

`treetime timetree` fails on `data/lassa/L/20` in the first backward pass, when it normalizes the combined time distribution of an internal node:

```
Error:
   0: When normalizing the time distribution of node 37
   1: Cannot normalize a distribution on [314.4969722119549, 1976.1488159980956]: its negative
      log-likelihood values contain NaN

Location: packages/treetime-distribution/src/distribution_core/distribution.rs:213
```

## Symptom and reproduction

```bash
treetime timetree --tree=data/lassa/L/20/tree.nwk --dates=data/lassa/L/20/metadata.tsv \
  --aln=data/lassa/L/20/aln.fasta.xz --output-all=<dir>
```

This is the `timetree/lassa/L/20/basic` case of `dev/smoke`.

The same abort occurs on `data/dengue/1000` (`timetree/dengue/1000/basic`, with `--name-column=genbank_accession`), in round 1 at node 240, on a distribution over [1846.41, 1966.49]. The `rust` branch fails on this dataset with `Cannot divide by empty distribution` instead, the error of [H-timetree-empty-distribution-division.md](H-timetree-empty-distribution-division.md).

Whether the dates of this dataset also have disjoint support is unknown until the `NaN` is removed.

### Simulated data

The timetree suite of [validation/README.md](../../validation/README.md) fails with this error in 13 runs of the core set, all with `--coalescent-opt` or `--coalescent-skyline`. Two of them, both `keep-skyline-fixed` on replicate `const/rep01` at alignment length $L = 10^4$:

- substitution rate 0.0003: the error chain starts with `When running round 2`
- substitution rate 0.001: the error chain has no round or node context, so it comes from a code path that adds none

Reproduce from `validation/`:

```bash
pixi run python simulate.py --coal const --rep 1
pixi run python run.py --suite timetree --set core --filter 'keep-skyline-fixed/const/rep01/L10000_'
```

Other failures of this suite produce different errors on other paths: [M-timetree-weak-signal-errors-differ-by-code-path.md](M-timetree-weak-signal-errors-differ-by-code-path.md).

## Impact and scope

- The command produces no output for this dataset.
- The case is excluded from every before/after output comparison, leaving this path without regression coverage.

## Mechanism

The combined backward message at the node is the product of the child messages, the coalescent prior, and the date constraint, so the `NaN` comes from one of those factors or from their product. In negative-log space a product is a sum, and `+inf + -inf` is `NaN`; a message with `-inf` or `NaN` values would produce this result. The source has not been traced.

When the `NaN` is removed, two contracts still meet here. The forward pass tolerates an empty posterior: when the message from the rest of the tree and the node's own date constraint have disjoint support, it warns and leaves the node undated. `fn collect_tree_events()` [packages/treetime/src/coalescent/events.rs#L72-L77](../../packages/treetime/src/coalescent/events.rs#L72-L77) then requires a time on every node to order the coalescent events, and errors on the first undated one. Its error text blames a stale coalescent model, which is misleading when the posterior is empty.

Deciding the intended behavior for that case is part of the fix: either an empty posterior must be resolved to a time (falling back to the constraint or to the parent's implied time), or the coalescent must define what an undated node contributes to the lineage count, or the disjoint-support condition must fail earlier with a diagnosis naming the conflicting dates.

## Related issues

- [H-timetree-empty-distribution-division.md](H-timetree-empty-distribution-division.md)
- [M-timetree-weak-signal-errors-differ-by-code-path.md](M-timetree-weak-signal-errors-differ-by-code-path.md)
- [H-timetree-mass-sizing-node-times-break-downstream-invariants.md](H-timetree-mass-sizing-node-times-break-downstream-invariants.md)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md)
