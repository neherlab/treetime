# Implement stochastic polytomy resolution

v1 implements only greedy deterministic polytomy resolution. v0 supports stochastic resolution via coalescent simulation for uncertainty quantification.

## Background

### Polytomies in phylogenetics

A polytomy (multifurcation) is a node with more than two children in a phylogenetic tree. Polytomies arise from:

- Soft polytomies: Insufficient data to resolve true bifurcating topology
- Hard polytomies: Genuine simultaneous divergence events (rare)

Most phylogenetic methods assume binary trees. Timetree inference requires resolving polytomies to estimate divergence times on internal branches.

### Resolution approaches

Greedy deterministic (implemented in v1):

1. Sort children by estimated divergence time
2. Iteratively merge closest pair using Brent optimization for branch lengths
3. Produces single resolved topology

Stochastic coalescent (v0 only):

1. Model polytomy as simultaneous sampling from a coalescent process
2. Sample binary resolution using Kingman coalescent
3. Repeat to generate distribution of topologies
4. Select best by likelihood or retain ensemble for uncertainty

### Kingman coalescent

The Kingman coalescent <a id="cite-1"></a>[Kingman 1982](<https://doi.org/10.1016/0304-4149(82)90011-4>) [[1](#ref-1)] describes the genealogical process of a sample
from a population. For k lineages, coalescent events occur at rate:

```
λ_k = k(k-1) / (2 * N_e)
```

where N_e is effective population size. The waiting time to next coalescent is exponentially distributed with rate λ_k. At each event, two lineages merge uniformly at random.

For polytomy resolution, the coalescent provides a principled model for how the unresolved lineages relate, given the time constraint from parent to children.

## v0 implementation

`generate_subtree()` (`#generate_subtree`) in [`packages/legacy/treetime/treetime/treetime.py#L872-L1011`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011):

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

- [`packages/treetime/src/timetree/optimization/polytomy.rs`](../../packages/treetime/src/timetree/optimization/polytomy.rs)
- Deterministic output
- No stochastic sampling
- No topology uncertainty quantification

## Impact

- Users needing topology uncertainty from polytomy resolution must use v0
- Greedy resolution produces a single "best" topology without confidence measure
- For trees with many polytomies, the deterministic choice may bias downstream analysis
- The star tree paradox <a id="cite-3"></a>[Lewis, Holder, and Holsinger 2005](https://doi.org/10.1080/10635150590924208) [[3](#ref-3)] shows that unresolved polytomies can bias inference toward the star topology in Bayesian analysis; stochastic resolution mitigates this by exploring the space of binary resolutions

## Algorithm detail

See [Stochastic Polytomy Resolution](../algo/unimplemented.md#stochastic-polytomy-resolution) for the full v0 algorithm walkthrough.

The TreeTime pipeline <a id="cite-2"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[2](#ref-2)] describes the coalescent-based polytomy resolution as part of the iterative timetree refinement.

## Cross-references

- [Polytomy resolution design](../reports/iterative-tree-refinement/7-polytomy-resolution.md): Detailed design for polytomy resolution in the iteration loop.

## References

1. <a id="ref-1"></a> Kingman, John F. C. 1982. "The Coalescent." _Stochastic Processes and Their Applications_ 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4 [↩](#cite-1)
2. <a id="ref-2"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-2)
3. <a id="ref-3"></a> Lewis, Paul O., Mark T. Holder, and Kent E. Holsinger. 2005. "Polytomies and Bayesian Phylogenetic Inference." _Systematic Biology_ 54(2):241-253. https://doi.org/10.1080/10635150590924208 [↩](#cite-3)

## Related issues

- Source: [N-timetree-stochastic-polytomy-unimplemented.md](../issues/N-timetree-stochastic-polytomy-unimplemented.md) -- delete after full resolution
