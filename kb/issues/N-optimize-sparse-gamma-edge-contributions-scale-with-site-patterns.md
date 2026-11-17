# Sparse gamma edge contributions scale with the number of site patterns

Without rate variation, the sparse edge contribution of the branch-length step groups the sites that are not variable on an edge by their state pair, so its cost grows with the variable sites of the edge. With `--gamma-categories K`, each site has its own posterior over the K categories, which weights its terms, so the mixture contribution keeps every distinct site pattern on every edge (`fn sparse_mixture_contribution()`, [packages/treetime/src/partition/marginal/mixture/contribution.rs#L71](../../packages/treetime/src/partition/marginal/mixture/contribution.rs#L71)). Sites with identical category likelihoods share one row, which leaves 5943 patterns of 85115 listed sites on `data/mpox/clade-ii/100`.

On `data/mpox/clade-ii/500` (sparse) with K = 4, the branch-length step of one optimize iteration takes 47.5 s against 0.41 s without gamma, and it is about three quarters of the iteration. The whole run takes 25 times as long as without gamma ([kb/decisions/discrete-gamma-rate-variation.md](../decisions/discrete-gamma-rate-variation.md#cost)).

> [!IMPORTANT]
> **Decision required.** The exact mixture objective of an edge is a sum over sites of the log of a K-term sum whose weights differ per site, so an exact evaluation reads every pattern. The options are: keep the exact cost and document it; merge the patterns of a state pair whose category weights agree within a tolerance, with a proven bound on the error of the branch-length objective; or evaluate the branch-length step of a mixture on a reduced set of patterns and keep the exact objective for the stop rule. Each approximation changes the branch lengths and needs an oracle comparison against the exact objective (dense mode, IQ-TREE).

## Validation

- The branch-length step of mpox/clade-ii/500 with K = 4 costs at most about K times the step without gamma, or the decision file records the exact cost as intended
