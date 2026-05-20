# Dense-sparse log-likelihood divergence

Dense and sparse marginal reconstruction produce different log-likelihoods for ~2.5% of random gap-free GTR configurations. The property test [`test_prop_marginal_dense_sparse_gap_free_consistency`](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L11) (`#test_prop_marginal_dense_sparse_gap_free_consistency`) detects two distinct populations. Under investigation.

- Population 1 (~97.5%): agrees to 0-3 ULPs (floating-point rounding)
- Population 2 (~2.5%): relative differences up to 7.6e-6 (billions of ULPs)

The test uses `max_relative = 1e-5` to accommodate the worst observed case.

## Investigation leads

Candidates for root cause (from test file comments):

1. Extreme eigenvalue ratios in GTR rate matrix Q causing numerical divergence
   in matrix exponentiation between the two code paths
2. Near-degenerate exchangeability matrices where small perturbations change
   the Fitch variable/invariant classification
3. Fitch invariant/variable classification diverging from likelihood-based
   variability assessment (some positions classified differently by parsimony vs probability)

The bimodal distribution (clean separation between populations 1 and 2) suggests a discrete trigger rather than continuous numerical drift.

4. Sparse EPS demotion: `combine_messages()` demotes variable sites to fixed
   when posterior peak exceeds `1 - EPS` (`EPS = 1e-4` at `packages/treetime/src/partition/marginal_helpers.rs#L15`). This is an approximation not present in the dense path, which retains full per-position probability vectors. The demotion replaces a position-specific posterior with a shared per-character vector, introducing a small error proportional to the distance between the position posterior and the shared vector.

## Excluded cause

Partial-IUPAC sparse reference-state fallback is covered separately by [`test_marginal_dense_sparse_ambiguous_r_reference_state_consistency`](../../packages/treetime/src/ancestral/__tests__/test_marginal_consistency.rs#L317) and the `flu/h3n2/500` site-451 regression. The remaining issue in this file is not the fixed ambiguous-`R` exact-state bug.
