# The reported clock rate is not the rate the reported dates were inferred at

Each refinement round infers the node times and then steps the clock rate (`fn update_clock_model` in
[round.rs](../../packages/treetime/src/timetree/round.rs)). Without `--time-marginal` (always or only-final) and
without `--confidence`, `fn refine_final_times` in [pipeline.rs](../../packages/treetime/src/timetree/pipeline.rs)
keeps the last round's time inference. The output therefore reports the rate after the last step, while its dates
come from the rate before it. With the default two rounds, a change to the last rate step changes the reported rate
but no date.

Found in the full validation of the rate window (2880 runs, commit fe7037bd4 against 26c1c8533): the median rate
error at $L=10^3$, $\mu=3\cdot10^{-4}$ (`keep-opt-est`) fell from 0.505 to 0.281, and not one run changed its date
RMSE. The gate of [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window)
ran with `--time-marginal only-final`, whose final pass uses the final rate.

> [!IMPORTANT]
> **Decision required.** Either a final time pass at the final rate in every mode (one more backward and forward pass,
> as `only-final` does), or report the rate the dates were inferred at and drop the last step. The first makes the
> output consistent and lets the last step count; the second costs nothing but wastes the last step's probes.

## Validation

For every run of the validation suite, the clock model's rate and the rate of the final time inference agree, and
the date RMSE follows the rate step.
