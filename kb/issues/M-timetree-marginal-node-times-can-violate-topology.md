# Marginal timetree node times can violate topology

The marginal forward pass selects the peak of each node's posterior independently. It does not jointly constrain the selected peaks to satisfy `child_time >= parent_time` in calendar coordinates, so an inferred edge can have a negative time length even when each node individually maximizes its marginal posterior.

## Current and reference behavior

[`set_likely_time()`](../../packages/treetime/src/timetree/inference/forward_pass.rs) assigns each node's `likely_time()` without consulting its parent. Later edge-time assignment derives the duration from the two selected node times.

V0 [`ClockTree._ml_t_marginal()`](../../packages/legacy/treetime/treetime/clock_tree.py) explicitly states that marginal reconstruction can produce negative branch lengths, then assigns the independent marginal peak and computes `clock_length` from the parent and child peaks. The parent-plus-epsilon clamp proposed in commit `542ac860c7cfa4bab6764aee1d1b3810a09eb54f` therefore does not match v0 behavior.

## Why a post-hoc clamp is not neutral

Changing only the stored node time to `parent_time + epsilon` would leave three representations inconsistent:

- the node posterior still peaks at the original unconstrained time;
- confidence intervals and likelihood-derived state still describe that posterior;
- edge time length and output date would describe the clamped value.

The proposed fixed epsilon also has no derivation from the time coordinate, grid spacing, posterior resolution, or data precision.

## Decision required

No implementation ticket is ready until the project chooses the desired contract:

- preserve v0 parity and permit/report negative marginal edge lengths;
- perform constrained inference so node posteriors and selected times satisfy topology together; or
- approve an output projection and define how posteriors, confidence intervals, likelihoods, and edge state are recomputed or labeled after projection.

## Related issues

- [M-distribution-support-boundary-semantics-unresolved.md](M-distribution-support-boundary-semantics-unresolved.md) - tail semantics can change posterior peaks but do not guarantee topology
- [M-timetree-coalescent-branch-length-clamp.md](M-timetree-coalescent-branch-length-clamp.md) - coalescent likelihood separately clamps negative edge lengths
