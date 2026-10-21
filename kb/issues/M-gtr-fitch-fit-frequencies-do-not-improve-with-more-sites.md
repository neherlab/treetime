# The equilibrium frequencies of the Fitch GTR fit do not improve with more sites

With `--model infer`, the sequence commands fit the GTR once, from the Fitch counts ([kb/decisions/dense-fitch-gtr-inference.md](../decisions/dense-fitch-gtr-inference.md)). On simulated data with a known model, the error of the exchangeabilities falls when the alignment grows, but the error of the equilibrium frequencies does not fall detectably.

## Evidence

`test_sim_ancestral_gtr_error_falls_with_more_sites` [packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs](../../packages/treetime/src/ancestral/__tests__/test_sim_ancestral_gtr_error_reduction.rs): 200 replicates of coalescent trees of 40 tips (8 batches of 5 over 2 years, constant `Tc` of 5), a random GTR with `mu = 0.4` per year, exact branch lengths, dense marginal reconstruction with `GtrModelName::Infer`, on the first 1000 and on all 10000 sites. Log ratio of the error at 10000 and at 1000 sites (slow tier, 2026-10-21):

- **RMS error of `log w`**: mean -0.034, standard deviation 0.117, p(mean < 0) = 3.6e-5, below the alpha of 5e-4. The error falls
- **Largest error of `pi`**: mean -0.031, standard deviation 0.218, p(mean < 0) = 0.024, above alpha. No detectable fall

An earlier measurement on the same simulation found the rate of the single fit at 1, the truth, and `pi` within about 0.02 of the true frequencies.

The fit is `fn infer_gtr_fitch()` [packages/treetime/src/partition/fitch/gtr_inference.rs#L13-L24](../../packages/treetime/src/partition/fitch/gtr_inference.rs#L13-L24), which counts the Fitch transitions and dwell times and calls `fn infer_gtr_from_counts()` [packages/treetime/src/gtr/infer_gtr/common.rs#L14-L87](../../packages/treetime/src/gtr/infer_gtr/common.rs#L14-L87).

> [!IMPORTANT]
> **Investigation required.** Whether the error of `pi` has a bias that does not shrink with the number of sites (from the Fitch root states, from the dwell times, or from the pseudo-count), or whether the test needs more replicates to show the fall. Measure the mean signed error of each frequency against the true `pi` at 1000 and at 10000 sites, and compare the frequencies of the Fitch fit with the composition of the simulated root sequence.

## Validation

- `test_sim_ancestral_gtr_error_falls_with_more_sites` passes in the slow tier and loses its `#[ignore]`
