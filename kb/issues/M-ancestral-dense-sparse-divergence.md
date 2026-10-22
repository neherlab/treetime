# Dense-sparse log-likelihood divergence

Dense and sparse marginal reconstruction produce different log-likelihoods beyond rounding. How often depends on the input. Two measures describe it: the share of cases that differ beyond rounding, and the share above the relative tolerance `1e-5` of the property tests:

- Small random inputs (4 taxa, 10 gap-free sites, random GTR parameters): about 2.5% of the cases differ beyond rounding, about 97.5% agree to a few ULPs, and no case is above `1e-5` (measurement below)
- Simulated data (coalescent trees of 6 to 20 tips, a random GTR, 300 to 1500 sites, true GTR and exact branch lengths): 993 of 1000 cases (99.3%) differ beyond rounding (relative difference above 1e-12), and 9 of 1000 cases (0.9%) are above `1e-5`, with a median relative difference of 2.5e-7 and a largest one of 1.0e-4 (see Measurement on simulated data)

The property test `test_prop_marginal_dense_sparse_gap_free_consistency` [packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L14-L21](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L14-L21) detects two distinct populations on the small random inputs. The test is `#[ignore]`d, so CI does not run this comparison. The cause is unknown.

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

The property test `test_prop_marginal_dense_sparse_simulated_agree_on_likelihood` [packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs#L29](../../packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs#L29) runs both representations on data from the simulation generator (`arb_simulated_case()`: coalescent trees of 6 to 20 tips, a random GTR, 300 to 1500 sites), with the true GTR and the exact branch lengths. It asserts the relative tolerance `1e-5` above (`epsilon = 0.0`) over 64 cases, and it is `#[ignore]`d because of this issue. The states comparison of the same file is a separate test, which [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md) covers.

Measurement: 1000 cases of `arb_simulated_case()` from a deterministic proptest runner, with the relative difference `|dense - sparse| / max(|dense|, |sparse|)` of the log-likelihoods:

- 9 of 1000 cases (0.9%) are above `1e-5`. With 64 new random cases per run, about 44% of the runs contain one such case. The committed regression case below is also above `1e-5`, and proptest runs it before the new cases, so every run of the ignored test fails
- Median 2.5e-7, 90th percentile 6.5e-7, 99th percentile 8.4e-6, largest 1.0e-4
- 7 cases agree to rounding (below 1e-12, smallest 8.1e-16), and 3 more are below 1e-9. Agreement to rounding, as in population 1 of the small random inputs, therefore also occurs on simulated data, but it is rare. The other 993 cases differ beyond rounding, with relative differences from above 1e-12 up to 1.0e-4, and 5.0e-4 in the regression case below. This range reaches beyond the 7.6e-6 of population 2 on the small random inputs
- Largest case: dense -596.3121422628079, sparse -596.3741350729575, for `CaseSpec { seed: 1315641483007467270, sampling: [SampleBatch { tau: 0.0, count: 2 }, SampleBatch { tau: 0.5, count: 3 }, SampleBatch { tau: 1.0, count: 2 }, SampleBatch { tau: 1.5, count: 3 }, SampleBatch { tau: 2.0, count: 4 }], schedule: Exponential(ExponentialTc { tc0: 5.0, scale: 2.5 }), model: RandomGtr, mu: 0.003, length: 356 }`

A failing case of the routine test run, kept in the regression file [packages/treetime/proptest-regressions/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.txt](../../packages/treetime/proptest-regressions/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.txt), gives dense -1139.3036028844738 and sparse -1139.8773205007508 (relative difference 5.0e-4), with no differing most likely state:

```text
CaseSpec { seed: 16728650536634984346, sampling: [SampleBatch { tau: 0.0, count: 3 }, SampleBatch { tau: 0.6666666666666666, count: 3 }, SampleBatch { tau: 1.3333333333333333, count: 4 }, SampleBatch { tau: 2.0, count: 2 }], schedule: Exponential(ExponentialTc { tc0: 5.0, scale: 2.5 }), model: RandomGtr, mu: 0.003, length: 761 }
```

The generator does not cause the difference: the generator before its numerical fixes (commit `eb110bac3`) gives the same two log-likelihoods for this case up to the last digits (-1139.303602884453 and -1139.877320500729).

The calibration of the posteriors is equal for both representations on one simulated tree of 40 tips and 5000 sites: the mean of wrong most likely states minus their posterior error probability per site is -0.0209 for dense and -0.0210 for sparse, with p-values 0.490 and 0.487 against a mean of 0 ([packages/treetime/src/ancestral/__tests__/test_sim_ancestral_calibration.rs](../../packages/treetime/src/ancestral/__tests__/test_sim_ancestral_calibration.rs)). The divergence of the log-likelihoods therefore does not show in the calibration of the posteriors at this size

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
