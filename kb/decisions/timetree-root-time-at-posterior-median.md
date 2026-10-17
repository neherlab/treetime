# The root's time is committed at its posterior median

The forward pass commits each node at the peak of its posterior, except the root, which it commits
at the posterior median. The root's committed time is the value the refinement loop reads for the
clock model, the coalescent timescale and the lineage counts.

**Type**: Intentional divergence from v0, which commits every node at its posterior peak.

**v1 location**: `fn propagate_distributions_forward` in
[`forward_pass.rs`](../../packages/treetime/src/timetree/inference/forward_pass.rs), with
`fn posterior_median` in [`max_posterior_region.rs`](../../packages/treetime/src/timetree/max_posterior_region.rs).

## Why

The root's posterior is the broadest in the tree: only its two branches and the coalescent prior
constrain it. When it has two close peaks, the peak moves from one to the other between rounds, and
the reconstruction and the clock fit follow, which pushes it back. With the lineage counts rebuilt
once ([timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md#one-rebuild-then-frozen)),
`data/ebola/20` with `--coalescent-opt` flipped its root between 2013.96 and 2014.11 every round
while rate and $T_c$ had converged. The median moves continuously with the posterior.

## Results

`data/ebola/20` converges in 7 rounds with `--coalescent-opt` and 5 with `--coalescent 2.56`.
On simulated data with the rebuilt lineage counts (10 rounds, 30 sets per cell), the median halves
the largest node-time change after round 3 at $L=10^4$, $\mu=3\cdot10^{-4}$ (0.67 to 0.33 years)
and the worst root date error at $L=10^3$, $\mu=3\cdot10^{-4}$ (34 to 17 years). Elsewhere it
changes little.

## Scope

Only the root. Other nodes keep the peak: their posteriors are narrower, and a median would move
them off the maximum of a skewed posterior near a hard branch-length boundary.
