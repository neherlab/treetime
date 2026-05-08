# Sparse fixed-position scalar rate approximation

When per-site rate variation is active, the sparse marginal reconstruction path applies per-site rates only to variable positions. Fixed (invariant) positions continue using the scalar rate `mu`. This is a bounded approximation that v0 avoids by rejecting the combination entirely.

v0 raises `TypeError('TreeAnc: sequence compression and site specific gtr models are incompatible!')` at `packages/legacy/treetime/treetime/treeanc.py:186-187` when both `is_site_specific` and `compress` are true. v1 allows the combination with the approximation described here.

## Scientific context

Felsenstein's pruning algorithm computes likelihood as a product over independent sites. Each site `a` contributes a factor `P_a(t) = exp(Q * mu_a * t)` to the per-branch transition matrix. When rates differ by site, positions are no longer exchangeable: two sites with the same observed nucleotide but different rates produce different parent probability distributions.

The sparse representation groups invariant positions by observed state: all fixed-A sites share one propagated profile `P(t)[:, A]`, weighted by their count. With uniform rates, this grouping is exact because `P(t)[:, A]` is identical for all fixed-A positions. With per-site rates, each fixed-A position has a different `P(t, r_a)[:, A]`, but the sparse structure has no mechanism to represent this variation.

## Impact

For invariant positions (same nucleotide across all leaves):

- The ancestral state at every node is unambiguously the observed nucleotide. Per-site rates do not change this assignment.
- The rate affects the transition probability column `P(t, r)[:, A]`, which influences the parent's probability distribution. A fast-evolving invariant site (r=5) produces a wider parent distribution than a slow site (r=0.1).
- The likelihood contribution per position depends on the rate. Using scalar `mu` instead of position-specific `r_a` introduces an error proportional to `|r_a - 1| * t` (the deviation from the mean rate times the branch length).

The per-position error at a fixed site is proportional to `|r_a - 1| * t`. This local error propagates through the tree via message passing and affects ALL node profiles, including at variable positions. A dense-vs-sparse comparison on a 16bp test alignment with gamma rates showed ~13% profile divergence at variable positions. The magnitude on real datasets depends on alignment length, tree depth, and rate heterogeneity.

The approximation affects likelihood magnitude and marginal profiles at all positions (not just fixed positions), because the fixed-position error accumulates through tree traversal.

## Limitations

- Log-likelihood values from sparse mode with per-site rates are approximate, not exact. Dense mode produces exact likelihoods.
- Parent probability distributions at fixed positions are slightly biased toward the mean-rate profile. This bias is largest for sites with extreme rate multipliers on long branches.
- GTR inference from sparse mutation counts does not account for per-site rates (a separate issue tracked in `../issues/M-gtr-per-site-rate-variation.md`).

## Potential improvements

- Group fixed positions by `(state, rate_category)` instead of `state` alone when discrete gamma categories are used. This preserves exchangeability within each category and makes the computation exact for the discrete approximation.
- Disable fixed-site aggregation entirely for site-rate models, storing per-position profiles. This eliminates the approximation at the cost of losing the memory advantage of sparse compression.

## References

- Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359
- Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00178256
- v0 guard: `packages/legacy/treetime/treetime/treeanc.py:186-187`
