# Polytomy zero-branch penalty differs from v0

v1 uses a linear time penalty for zero-mutation branches introduced during
polytomy resolution. v0 uses `zero_branch_slope = gtr.mu * data.full_length`,
scaling the penalty by the expected number of mutations per unit time.

- v1: `PolytomyCostFn::cost()` (`#PolytomyCostFn`) at
  [`packages/treetime/src/commands/timetree/optimization/polytomy.rs#L308-L309`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L308-L309)
- v0: `zero_branch_slope` at
  [`packages/legacy/treetime/treetime/treetime.py#L726`](../../packages/legacy/treetime/treetime/treetime.py#L726),
  used in `_c_gain()` (`#_c_gain`) at
  [`packages/legacy/treetime/treetime/treetime.py#L744`](../../packages/legacy/treetime/treetime/treetime.py#L744)

## Background

Polytomy resolution introduces a new internal node between a parent and two of
its children. The new branch from parent to the internal node has zero observed
mutations. A penalty discourages placing the internal node too close to the
parent (trivially short branches) or too far (implausibly long zero-mutation
branches).

## v0: rate-scaled penalty

```python
zero_branch_slope = self.gtr.mu * self.data.full_length
cg_new = -zero_branch_slope * (parent.time_before_present - t)
```

The penalty scales linearly with time difference, weighted by `mu * L` (expected
substitutions per unit time for the full alignment). Longer alignments and
higher mutation rates produce stronger penalties, making the cost proportional
to the expected number of mutations on the zero-length branch.

## v1: unscaled linear penalty

```rust
let zero_branch_penalty = new_branch_to_parent;
```

The penalty is the bare time difference, without scaling by mutation rate or
sequence length. This makes polytomy resolution independent of the
substitution model parameters but underweights the penalty for long alignments
with high mutation rates, and overweights it for short alignments.

## Impact

Different penalties affect which child pairs are selected for merging and where
the internal node is placed. The effect is largest for datasets with many
polytomies and large `mu * L` values, where v0's rate-scaled penalty is
much larger than v1's unscaled penalty.
