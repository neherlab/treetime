# A coalescent time scale in substitutions per site gives a higher clock rate than the same time scale in years

With a fixed coalescent time scale and an estimated clock rate, `timetree` converts a time scale given in substitutions per site into years at the current clock rate in every round. Given the true $T_c$ in substitutions per site at the true rate, the run should settle at the clock rate of the same run with the true $T_c$ in years. After the divergence floor of the start rate ([kb/decisions/timetree-start-rate-falls-back-below-a-divergence-floor.md](../decisions/timetree-start-rate-falls-back-below-a-divergence-floor.md)), the rate of the run in substitutions per site is about 3% higher, and the simulation test `test_sim_timetree_tc_in_substitutions_per_site_settles_with_the_estimated_rate` ([packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs](../../packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs), `just test-slow-rs`) fails.

## Evidence

The test runs 633 simulated replicates of the standard timetree cell (`fn timetree_cell_params()`: `--keep-root`, `--clock-filter 0`, 2 rounds) and requires the mean log ratio of the two clock rates to lie within $\pm 0.03$ at a false-failure rate of $10^{-3}$. Result per commit:

- before the start-rate floor (with the shrinking probe window of the rate step): passes
- with the floor: fails, mean $0.0292$, standard deviation $0.126$, interval $[0.0136, 0.0448]$
- with the floor and the denser-end anchor of resampled functions: fails, mean $0.0307$, standard deviation $0.125$, interval $[0.0153, 0.0461]$; the later commits leave these numbers unchanged

The floor changes no output of the `timetree` smoke cases. It changes a run only where the regression slope is positive and at most the floor, so on these simulated cells it replaces the slope by the Gaussian rate in some replicates. The first clock fit does not read the coalescent time scale, so both runs of a replicate start from the same rate.

> [!IMPORTANT]
> **Investigation required.** Find why the two runs end at different rates after a start from the Gaussian rate: record per replicate whether the floor fired, the start rate, and the rate after each round of both runs, and repeat with more rounds to see whether the difference shrinks. A difference that shrinks with the rounds means that 2 rounds do not settle the conversion of the time scale; one that stays means that the two time scales define different fixed points.

## Validation

`just test-slow-rs -E 'test(/tc_in_substitutions_per_site_settles/)'` passes.
