# The equilibrium frequencies of the Fitch GTR fit do not improve with more sites

With `--model infer`, the sequence commands fit the GTR once, from the Fitch counts ([kb/decisions/dense-fitch-gtr-inference.md](../decisions/dense-fitch-gtr-inference.md)). On simulated data with a known model, the error of the exchangeabilities falls when the alignment grows, but the error of the equilibrium frequencies does not fall detectably.

## Evidence

`test_sim_ancestral_gtr_error_falls_with_more_sites` [packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs](../../packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs): 200 replicates of coalescent trees of 40 tips (8 batches of 5 over 2 years, constant `Tc` of 5), a random GTR with `mu = 0.4` per year, exact branch lengths, dense marginal reconstruction with `GtrModelName::Infer`, on the first 1000 and on all 10000 sites. Log ratio of the error at 10000 and at 1000 sites (slow tier, 2026-10-21):

- **RMS error of `log w`**: mean -0.034, standard deviation 0.117, p(mean < 0) = 3.6e-5, below the alpha of 5e-4. The error falls
- **Largest error of `pi`**: mean -0.031, standard deviation 0.218, p(mean < 0) = 0.024, above alpha. No detectable fall

An earlier measurement on the same simulation found the rate of the single fit at 1, the truth, and `pi` within about 0.02 of the true frequencies.

The fit is `fn infer_gtr_fitch()` [packages/treetime/src/partition/fitch/gtr_inference.rs#L12-L23](../../packages/treetime/src/partition/fitch/gtr_inference.rs#L12-L23), which counts the Fitch transitions and dwell times and calls `fn infer_gtr_from_counts()` [packages/treetime/src/gtr/infer_gtr/common.rs#L16-L25](../../packages/treetime/src/gtr/infer_gtr/common.rs#L16-L25).

> [!IMPORTANT]
> **Investigation required.** Whether the error of `pi` has a bias that does not shrink with the number of sites (from the Fitch root states, from the dwell times, or from the pseudo-count), or whether the test needs more replicates to show the fall. Measure the mean signed error of each frequency against the true `pi` at 1000 and at 10000 sites, and compare the frequencies of the Fitch fit with the composition of the simulated root sequence.

## Validation

- `test_sim_ancestral_gtr_error_falls_with_more_sites` passes in the slow tier and loses its `#[ignore]`

## Proposed implementation plan

The test exists and carries `#[ignore = "kb/issues/M-gtr-fitch-fit-frequencies-do-not-improve-with-more-sites.md"]`; `fn infer_gtr_fitch()` and `fn infer_gtr_from_counts()` are at the cited lines. The measured numbers come from earlier runs and were not measured again on the current code.

1. Investigation first (the block above): add a temporary diagnostic to the simulation test that writes, per replicate, the signed error of each `pi` component at 1000 and 10000 sites, the composition of the simulated root sequence, and the Fitch root-state counts (`counts.root_state` of `fn get_mutation_counts_fitch()`). Run it with `just test-slow-rs test_sim_ancestral_gtr_error_falls_with_more_sites --run-ignored all` (slow recipe: ask the user first).
2. Decide by the result:
   - **Bias that does not shrink** (mean signed error of `pi` stays away from 0 at 10000 sites): trace its source in `fn get_mutation_counts_fitch()` (root states, dwell times on floored branch lengths, see [N-gtr-transition-counting-uses-floored-branch-lengths.md](N-gtr-transition-counting-uses-floored-branch-lengths.md)) and in the pseudo-count of `InferGtrOptions::default()`; fix the cause in the counting code. Any change of the estimator is a scientific choice and needs user approval and an update of [kb/decisions/dense-fitch-gtr-inference.md](../decisions/dense-fitch-gtr-inference.md).
   - **Noise only** (no bias, large variance): raise the replicate count or the site ratio of the test until the power is enough at alpha 5e-4, and keep the oracle (the true simulated `pi`).
3. Remove the `#[ignore]` when the test passes, and delete this issue.

Tests: the existing `test_sim_ancestral_gtr_error_falls_with_more_sites` (oracle: the true GTR of `treetime-simulation`). Add a unit test of the counting fix, if any, on a small hand-made tree with a hand-computed `pi` estimate.

Validation: `just test-rs`, then the slow-tier test on user request.
