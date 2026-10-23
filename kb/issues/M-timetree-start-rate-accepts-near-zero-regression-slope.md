# Timetree start rate accepts a near-zero positive regression slope

With an estimated clock rate, timetree starts from the slope of the root-to-tip regression whenever that slope is positive, and from the rate of the Gaussian timetree only when it is not ([packages/treetime/src/clock/clock_regression.rs#L151-L162](../../packages/treetime/src/clock/clock_regression.rs#L151-L162), called by `fn estimate_initial_clock()` [packages/treetime/src/timetree/pipeline.rs#L278-L304](../../packages/treetime/src/timetree/pipeline.rs#L278-L304) and in every round [packages/treetime/src/timetree/round.rs#L307](../../packages/treetime/src/timetree/round.rs#L307)). On a tree without root-to-tip signal the slope scatters around zero, so a small change in the input flips between a sensible start (the Gaussian rate) and a start rate orders of magnitude too low. The refinement loop cannot recover from the low start: the first round moves the rate by at most 2-fold ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window)). The node dates then reach thousands of years before the samples, and the run either exceeds `--max-grid-points` or panics in the node-data writer ([H-timetree-negative-numdate-panics-in-date-string-conversion.md](H-timetree-negative-numdate-panics-in-date-string-conversion.md)).

## Reproduction

Messy data of the validation suite (`validation/degrade.py`; branch `messy-validation2`, commit `0ddc02058` with the fixes up to `046ff2977`), set `const/rep03/L1000_mu0.003` of `data_messy/` (400 tips, true rate $3 \cdot 10^{-3}$), from `validation/`:

```bash
../.out/treetime timetree -t data_messy/const/rep03/L1000_mu0.003/inputs/realized.nwk \
  -a data_messy/const/rep03/L1000_mu0.003/tips_observed.fasta.xz \
  -m data_messy/const/rep03/tip_dates_observed.tsv --date-column date --clock-filter 0 --keep-root --max-iter 0 \
  --output-clock-model clock.json
```

fails with `Time inference needs a grid of 1246429 points, more than the limit of 1000000 (--max-grid-points). The estimated clock rate 2.29e-5 may be far from the true rate`. With `--coalescent opt` or `--coalescent skyline` and the default iterations, the run completes the time inference and panics in `fn refine_dates()` with `cannot convert float seconds to Duration: value is negative`.

The tree is the same in the clean and the messy set. Crossing the inputs on this set, `--keep-root`, estimated rate:

| alignment | tip dates | regression slope (`treetime clock --keep-root`) | timetree start rate (`--max-iter 0`) | timetree final rate |
| --- | --- | --- | --- | --- |
| clean | true | $-4.70 \cdot 10^{-5}$ | $3.95 \cdot 10^{-3}$ (Gaussian fallback) | $2.97 \cdot 10^{-3}$ |
| messy | true | | $3.95 \cdot 10^{-3}$ | $2.88 \cdot 10^{-3}$ |
| clean | observed | $+2.29 \cdot 10^{-5}$ | grid limit | grid limit |
| messy | observed | | grid limit | grid limit |

The observed dates (3% of tips undated, 5% dated to the month only) turn a slightly negative slope into a slightly positive one; the alignment does not matter. The failure happens before the refinement loop, so the in-loop regression of undated tips ([M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md](M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md)) is not its cause. With `--reroot least-squares` the same set gives $3.02 \cdot 10^{-3}$. At $\mu = 10^{-3}$ the same replicate fails the same way (start rate $3.2 \cdot 10^{-6}$).

## Frequency

Timetree `core` set of the suite, 10 tips per sample time, 240 sets per flavor, the six variants with estimated rate (1440 runs per flavor):

- messy data: 4 failed runs, all `keep-*-est` on `const/rep03` at $L = 10^3$ (2 grid-limit errors, 2 panics)
- rate below half the true rate in completed runs: 22 on clean data, 19 on messy data, most at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ (about one substitution per genome over the four years of sampling), where weak signal is expected. `const/rep03` is the exception: on clean data `keep-none-est` underestimates by 78% at $\mu = 10^{-3}$ and by 51% at $\mu = 10^{-2}$

## v0

`ClockTree.get_clock_model()` uses the regression slope as it is and fails only on a non-finite slope ([packages/legacy/treetime/treetime/clock_tree.py#L295-L311](../../packages/legacy/treetime/treetime/clock_tree.py#L295-L311)); it has no fallback for a non-positive slope. The Gaussian fallback of v1 is new. v0 has not been run on this set.

## Prior measurement

[kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window) tried the Gaussian-timetree start rate in place of the regression for all runs: it fixed `const/rep03` (then at $\mu = 10^{-2}$, a start at half the true rate) but raised the date RMSE at $L = 10^3$, $\mu = 3 \cdot 10^{-3}$ from 0.155 to 0.191 years, so the regression start was kept. That measurement predates the messy data and the hard failures above.

> [!IMPORTANT]
> **Decision required.** The start rate switches at a slope of exactly zero, although a slope near zero carries no more information than a slightly negative one. Options:
>
> - **Gaussian start below a signal threshold:** use the Gaussian-timetree rate when the regression is weak, for example when the slope is within its standard error of zero or the explained variance is below a bound, and the regression slope otherwise. Keeps the regression start, and its RMSE advantage, where the signal is clear
> - **Gaussian start always:** as measured in the decision above, at the cost recorded there
> - **Larger first-round rate window:** keep the start and let the first round move the rate by more than 2-fold when the probes keep rising, so a start too low by orders of magnitude recovers. Measured on the reproduction with `keep-opt-est`: a 10-fold cap, the full step and a cube-root parabola end at $4.6 \cdot 10^{-5}$, because $T_c$ fits the stretched tree and $\log Z$ under that prior peaks near the current rate ([kb/reports/timetree-rate-step-coordinate-and-guards.md](../reports/timetree-rate-step-coordinate-and-guards.md#far-off-start))
> - **Reject the input:** return an input error that names the weak date signal when the start rate is implausible, as `NoPositiveRateRoot` does for the regression ([M-timetree-weak-signal-errors-differ-by-code-path.md](M-timetree-weak-signal-errors-differ-by-code-path.md))
>
> Evidence: the table and frequencies above; the RMSE cost of the Gaussian start in the decision above.

## Validation

The timetree suite on `data_messy/` completes `keep-*-est` on `const/rep03` with a rate within the spread of the other replicates, and the date RMSE of the clean `core` set at $L = 10^3$, $\mu = 3 \cdot 10^{-3}$ does not rise above that of the regression start.
