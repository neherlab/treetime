# Stochastic polytomy resolution not implemented

v1 implements only greedy deterministic polytomy resolution. v0 supports stochastic
resolution via coalescent simulation for uncertainty quantification.

## Background

### Polytomies in phylogenetics

A polytomy (multifurcation) is a node with more than two children in a phylogenetic
tree. Polytomies arise from:

- **Soft polytomies**: Insufficient data to resolve true bifurcating topology
- **Hard polytomies**: Genuine simultaneous divergence events (rare)

Most phylogenetic methods assume binary trees. Timetree inference requires resolving
polytomies to estimate divergence times on internal branches.

### Resolution approaches

**Greedy deterministic** (implemented in v1):

1. Sort children by estimated divergence time
2. Iteratively merge closest pair using Brent optimization for branch lengths
3. Produces single resolved topology

**Stochastic coalescent** (v0 only):

1. Model polytomy as simultaneous sampling from a coalescent process
2. Sample binary resolution using Kingman coalescent
3. Repeat to generate distribution of topologies
4. Select best by likelihood or retain ensemble for uncertainty

### Kingman coalescent

The Kingman coalescent (Kingman 1982) describes the genealogical process of a sample
from a population. For k lineages, coalescent events occur at rate:

```
λ_k = k(k-1) / (2 * N_e)
```

where N_e is effective population size. The waiting time to next coalescent is
exponentially distributed with rate λ_k. At each event, two lineages merge uniformly
at random.

For polytomy resolution, the coalescent provides a principled model for how the
unresolved lineages relate, given the time constraint from parent to children.

## v0 implementation

`generate_subtree()` (`#generate_subtree`) in
[`packages/legacy/treetime/treetime/treetime.py#L872-L1011`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011):

1. Sort polytomy children by `time_before_present`
2. Initialize coalescent rate from `merger_model` (if available) or dummy rate
3. Loop until binary or time exhausted:
   - Branches without pending mutations are "ready to coalesce"
   - Sample waiting time from exponential distribution
   - Either remove a mutation from a branch or coalesce two ready branches
   - Coalescence picks two lineages uniformly at random
4. Remaining uncoalesced branches become direct children of parent

Key parameters:

- `merger_model.branch_merger_rate(t)`: Time-varying coalescent rate from skyline
- `mutation_rate = gtr.mu * L`: Substitution rate affects event timing
- Uses TreeTime's RNG (`self.rng.exponential`, `self.rng.choice`)

## v1 status

v1 uses greedy pairwise merging with Brent optimization:

- [`packages/treetime/src/commands/timetree/optimization/polytomy.rs`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs)
- Deterministic output
- No stochastic sampling
- No topology uncertainty quantification

## Impact

- Users needing topology uncertainty from polytomy resolution must use v0
- Greedy resolution produces a single "best" topology without confidence measure
- For trees with many polytomies, the deterministic choice may bias downstream analysis

## Algorithm detail

See [Stochastic Polytomy Resolution](../port-algo-inventory/unimplemented.md#stochastic-polytomy-resolution) for the full v0 algorithm walkthrough.

## References

- Kingman, J.F.C. (1982). "The coalescent." Stochastic Processes and their
  Applications, 13(3):235-248.
- Sagulenko, P., Puller, V., & Neher, R.A. (2018). "TreeTime: Maximum-likelihood
  phylodynamic analysis." Virus Evolution, 4(1):vex042.
