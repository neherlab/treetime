# Dense-sparse log-likelihood divergence

Dense and sparse marginal reconstruction produce different log-likelihoods for ~2.5% of random gap-free GTR configurations. The property test `test_prop_marginal_dense_sparse_gap_free_consistency` [packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L14-L21](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L14-L21) detects two distinct populations. The test is `#[ignore]`d, so CI does not run this comparison. The cause is unknown.

- Population 1 (~97.5%): agrees to 0-3 ULPs, `rel_diff < 5e-16` (floating-point rounding from a different order of operations)
- Population 2 (~2.5%): relative differences up to 7.6e-6 (billions of ULPs); these are algorithmic differences, not rounding

Measurement: one run of 200 proptest cases (the test runs 30), 4 taxa, 10-site gap-free alignments, random GTR parameters. Dense and sparse marginal reconstruction ran on identical inputs, and each case recorded the absolute difference, the relative difference, and the ULP distance. Observed cases of population 2:

| abs_diff | rel_diff | ULPs  |
| -------- | -------- | ----- |
| 3.0e-13  | 3.4e-15  | 21    |
| 3.3e-6   | 3.4e-8   | 234M  |
| 1.8e-5   | 1.8e-7   | 1.3B  |
| 3.0e-5   | 4.7e-7   | 4.2B  |
| 4.7e-5   | 5.8e-7   | 3.3B  |
| 5.6e-4   | 7.6e-6   | 39.5B |

The test uses `max_relative = 1e-5` to accommodate the worst observed case (7.6e-6). A ULP tolerance does not fit the bimodal distribution. `epsilon = 0.0` turns off the absolute tolerance of `relative_eq!`, so the comparison is purely relative.

## Measurement on simulated data

`test_prop_marginal_dense_sparse_simulated_agree_on_states_and_likelihood` [packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs](../../packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs) runs both representations on data from the simulation generator (coalescent trees of 6 to 20 tips, a random GTR, 300 to 1500 sites), with the true GTR and the exact branch lengths. It asserts the relative tolerance `1e-5` above. One run of its 64 cases measured:

- Every case differs by more than rounding: relative differences of the log-likelihoods from 1.2e-9 to 5.3e-6. On realistic data the two populations above do not separate; all cases fall into population 2
- The most likely states of the internal nodes agree at all 564950 compared node sites (dense posterior gap above `1e-3`), and at every excluded near-tie site
- The calibration of the posteriors is equal for both representations on one simulated tree of 40 tips and 5000 sites: the mean of wrong most likely states minus their posterior error probability per site is -0.0209 for dense and -0.0210 for sparse, with p-values 0.490 and 0.487 against a mean of 0 ([packages/treetime/src/ancestral/__tests__/test_sim_ancestral_calibration.rs](../../packages/treetime/src/ancestral/__tests__/test_sim_ancestral_calibration.rs)). The divergence of the log-likelihoods therefore does not show in the calibration of the posteriors at this size

## Investigation leads

Candidates for root cause (from test file comments):

1. Extreme eigenvalue ratios in GTR rate matrix Q causing numerical divergence
   in matrix exponentiation between the two code paths
2. Near-degenerate exchangeability matrices where small perturbations change
   the Fitch variable/invariant classification
3. Fitch invariant/variable classification diverging from likelihood-based
   variability assessment (some positions classified differently by parsimony vs probability).
   A position where all tips share one state is invariant to Fitch parsimony, but its
   probability vector under the GTR model still depends on the transition probabilities
   of each branch length

The bimodal distribution (clean separation between populations 1 and 2) suggests a discrete trigger rather than continuous numerical drift.

To find the trigger, compare intermediate results of the two paths on a population 2 case: per-node partial likelihoods, and the propagation along each branch through `P(t) = exp(Qt)`. If the cause is a correctness defect, fix it and tighten the test tolerance to the ULP level. If the difference is inherent to the sparse approximation, document its error bound.

4. Sparse EPS demotion: `fn combine_messages()` [packages/treetime/src/partition/marginal/sparse/message.rs#L22](../../packages/treetime/src/partition/marginal/sparse/message.rs#L22) demotes variable sites to fixed
   when the posterior peak exceeds `1 - EPS` (`const EPS: f64 = 1e-4` [packages/treetime/src/partition/marginal/sparse/message.rs#L20](../../packages/treetime/src/partition/marginal/sparse/message.rs#L20), applied through `is_site_resolved()` at [packages/treetime/src/partition/marginal/sparse/message.rs#L84](../../packages/treetime/src/partition/marginal/sparse/message.rs#L84)). This is an approximation not present in the dense path, which retains full per-position probability vectors. The demotion replaces a position-specific posterior with a shared per-character vector, introducing a small error proportional to the distance between the position posterior and the shared vector.

   The same demotion also drove a reconstruction defect, now fixed and recorded in [kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md). That defect affected emitted sequences only and left log-likelihoods untouched, so it neither explains nor eliminates this lead - the demotion's effect on log-likelihood still needs evaluating on its own.

## Excluded cause

Partial-IUPAC sparse reference-state fallback is covered separately by [`test_marginal_dense_sparse_ambiguous_r_reference_state_consistency`](../../packages/treetime/src/ancestral/__tests__/test_marginal_consistency.rs#L257) and the `flu/h3n2/500` site-451 regression. The remaining issue in this file is not the fixed ambiguous-`R` exact-state bug.
