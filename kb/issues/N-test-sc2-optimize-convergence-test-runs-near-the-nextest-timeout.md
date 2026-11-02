# The sc2 optimize convergence test runs near the nextest timeout

`test_convergence_sc2_sparse_converges_on_sc2_2844` ([packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L22-L76](../../packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L22-L76)) runs the sparse optimize loop on `sc2/2844` (2504 tips) for up to 50 iterations and asserts that it stops before the limit. Commit `8d32da409` added it as the real-dataset guard of `0d75d2636`, which made the convergence loop robust to the two-cycle of the sparse objective. `just test-rs` builds it in the `dev` profile, in which the workspace crates are not optimized.

The default nextest profile terminates a test after five periods of 60 s, 300 s in total (`slow-timeout` in `.config/nextest.toml`). In the full `just test-rs` runs of one orchestration on the shared build host, the test took 223 to 248 s. At a load average of 16 to 18 it exceeded 300 s once and failed with a timeout; run alone at a load average of 8.7, it passed in 227.6 s. The margin of about 20% makes the routine suite fail on a loaded host without a defect.

> [!IMPORTANT]
> **Decision required.** Three options keep the guard and remove the timeouts:
>
> - **Slow tier**: run the test with `just test-slow-rs` (release build, 40-minute stop). The `slow` profile selects tests by a `test_sim_` path segment, so the test needs that prefix or the filter of both profiles needs a second pattern. No CI job runs `test-slow-rs`, so CI would no longer run the guard
> - **Per-test timeout**: a `[[profile.default.overrides]]` entry in `.config/nextest.toml` that gives this test a longer `slow-timeout`. The test stays in the routine suite and in CI, and each `test-rs` run keeps paying its 4 minutes
> - **Smaller dataset**: a subset of `sc2/2844`, or another dataset, on which the loop without the convergence conditions of `0d75d2636` still enters the two-cycle. Fast in the routine suite, but the property must be shown again on the new input

> [!IMPORTANT]
> **Investigation required.** For the third option, find a dataset that reproduces the two-cycle when the convergence conditions of `0d75d2636` are turned off, and measure the run time of the test on it in the `dev` profile.
