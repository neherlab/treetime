# The timetree start rate is too high on deep trees with a short sampling window

With an estimated clock rate, `timetree` starts from the rate of the Gaussian timetree when the slope of the root-to-tip regression is at most the divergence floor $0.1\,\bar d/(t_\text{max} - t_\text{min})$ ([kb/decisions/timetree-start-rate-falls-back-below-a-divergence-floor.md](../decisions/timetree-start-rate-falls-back-below-a-divergence-floor.md)). Since $\bar d/(t_\text{max} - t_\text{min}) \approx \mu \cdot \text{depth}/\text{span}$, the floor lies above the true rate when the root is more than about ten sampling spans deep, and then fires on every input. The Gaussian substitution rate that it falls back to is biased high on such trees: its conversion $\mu = \eta D/(D+M+N)$ assumes one hazard for the branching and sampling events of all branches, while a dense recent sample has a much higher sampling hazard on terminal branches than the deep internal branches have (Section 15.4 of the maintainers' derivation of the Gaussian timetree).

## Evidence

- Keeping only the last 0.3 years of samples of 20 trees of `data_n30` (120 tips, root 18 to 100 sampling spans deep): the floor lay at 1.8 to 10 times the true rate and fired on every set, although the regression slope was $-0.7$ to $2.6$ times the true rate. The Gaussian rate was 1.4 to 2.6 times the true rate. Without a coalescent prior the final rates were unchanged; with `--coalescent opt` the median $\lvert\ln\rvert$ final rate error rose from 0.115 to 0.285, with no failures
- On the trees of the validation suite the floor is 0.2 (median) to 0.86 of the true rate, so it fires there only on weak signal

## Fix direction

A two-class hazard in the conversion $\mu = \eta D/(D+M+N)$, one for terminal and one for internal branches, is the remedy that the derivation documents. It is not implemented.

## Validation

On the shortened `data_n30` sets, the start rate lies within the spread of the regression slope around the true rate, and the median $\lvert\ln\rvert$ final rate error with `--coalescent opt` is not above that without the floor (0.115).

## Before the divergence floor

The measurements of this section ran on code before the divergence floor (`origin/rust` at `9b10ffd3d`). There, `timetree` started from the regression slope whenever it was positive, and from the rate of the Gaussian timetree only when it was not. On a tree without root-to-tip signal the slope scatters around zero, so a small change in the input flipped between a sensible start (the Gaussian rate) and a start rate orders of magnitude too low. The refinement loop cannot recover from the low start: the first round moves the rate by at most 2-fold ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window)). The node dates then reached thousands of years before the samples, and the run either exceeded `--max-grid-points` or panicked in the node-data writer ([H-io-negative-numdate-panics-in-date-string-conversion.md](H-io-negative-numdate-panics-in-date-string-conversion.md)). The divergence floor fixes this case: the start there is the Gaussian rate, and the runs end at 0.96 to 1.01 of the true rate (see the decision).

### Reproduction

Messy data of the validation suite (`validation/degrade.py`), set `const/rep03/L1000_mu0.003` of `data_messy/` (400 tips, true rate $3 \cdot 10^{-3}$), from `validation/`:

```bash
../.out/treetime timetree -t data_messy/const/rep03/L1000_mu0.003/inputs/realized.nwk \
  -a data_messy/const/rep03/L1000_mu0.003/tips_observed.fasta.xz \
  -m data_messy/const/rep03/tip_dates_observed.tsv --date-column date --clock-filter 0 --keep-root --max-iter 0 \
  --output-clock-model clock.json
```

failed with `Time inference needs a grid of 1246429 points, more than the limit of 1000000 (--max-grid-points). The estimated clock rate 2.29e-5 may be far from the true rate`. With `--coalescent opt` or `--coalescent skyline` and the default iterations, the run completed the time inference and panicked in `fn refine_dates()` with `cannot convert float seconds to Duration: value is negative`.

The tree is the same in the clean and the messy set. Crossing the inputs on this set, `--keep-root`, estimated rate:

| alignment | tip dates | regression slope (`treetime clock --keep-root`) | timetree start rate (`--max-iter 0`) | timetree final rate |
| --- | --- | --- | --- | --- |
| clean | true | $-4.70 \cdot 10^{-5}$ | $3.95 \cdot 10^{-3}$ (Gaussian fallback) | $2.97 \cdot 10^{-3}$ |
| messy | true | | $3.95 \cdot 10^{-3}$ | $2.88 \cdot 10^{-3}$ |
| clean | observed | $+2.29 \cdot 10^{-5}$ | grid limit | grid limit |
| messy | observed | | grid limit | grid limit |

The observed dates (3% of tips undated, 5% dated to the month only) turn a slightly negative slope into a slightly positive one; the alignment does not matter. The failure happens before the refinement loop, so the in-loop regression of undated tips ([M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md](M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md)) is not its cause. With `--reroot least-squares` the same set gives $3.02 \cdot 10^{-3}$. At $\mu = 10^{-3}$ the same replicate fails the same way (start rate $3.2 \cdot 10^{-6}$).

### Frequency

Timetree `core` set of the suite, 10 tips per sample time, 240 sets per flavor, the six variants with estimated rate (1440 runs per flavor):

- messy data: 4 failed runs, all `keep-*-est` on `const/rep03` at $L = 10^3$ (2 grid-limit errors, 2 panics)
- 30 tips per sample time, messy data: 3 failed runs, `keep-{none,opt,skyline}-est` on `const/rep05` at $L = 10^3$, $\mu = 10^{-3}$, all from the start rate $1.24 \cdot 10^{-6}$ (grid limit). At 100 tips per sample time the failures have another cause ([M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md](M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md))
- rate below half the true rate in completed runs: 22 on clean data, 19 on messy data, most at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ (about one substitution per genome over the four years of sampling), where weak signal is expected. `const/rep03` is the exception: on clean data `keep-none-est` underestimates by 78% at $\mu = 10^{-3}$ and by 51% at $\mu = 10^{-2}$

### v0

`ClockTree.get_clock_model()` uses the regression slope as it is and fails only on a non-finite slope ([packages/legacy/treetime/treetime/clock_tree.py#L295-L311](../../packages/legacy/treetime/treetime/clock_tree.py#L295-L311)); it has no fallback for a non-positive slope. The Gaussian fallback of v1 is new. v0 has not been run on this set.

### Prior measurement

[kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window) tried the Gaussian-timetree start rate in place of the regression for all runs: it fixed `const/rep03` (then at $\mu = 10^{-2}$, a start at half the true rate) but raised the date RMSE at $L = 10^3$, $\mu = 3 \cdot 10^{-3}$ from 0.155 to 0.191 years, so the regression start was kept. That measurement predates the messy data and the hard failures above.

### Options considered

The divergence floor is the first of these options; the others stay as measured alternatives:

- **Gaussian start below a signal threshold:** use the Gaussian-timetree rate when the regression is weak, for example when the slope is within its standard error of zero or the explained variance is below a bound, and the regression slope otherwise. Keeps the regression start, and its RMSE advantage, where the signal is clear
- **Gaussian start always:** as measured in the decision of the prior measurement above, at the cost recorded there
- **Larger first-round rate window:** keep the start and let the first round move the rate by more than 2-fold when the probes keep rising, so a start too low by orders of magnitude recovers. Measured on the reproduction with `keep-opt-est`: a 10-fold cap, the full step and a cube-root parabola end at $4.6 \cdot 10^{-5}$, because $T_c$ fits the stretched tree and $\log Z$ under that prior peaks near the current rate ([kb/reports/timetree-rate-step-coordinate-and-guards.md](../reports/timetree-rate-step-coordinate-and-guards.md#far-off-start))
- **Reject the input:** return an input error that names the weak date signal when the start rate is implausible, as `NoPositiveRateRoot` does for the regression ([M-timetree-weak-signal-errors-differ-by-code-path.md](M-timetree-weak-signal-errors-differ-by-code-path.md))

### Criterion of this case

The run on `data_messy/` was to complete `keep-*-est` on `const/rep03` with a rate within the spread of the other replicates, and the date RMSE of the clean `core` set at $L = 10^3$, $\mu = 3 \cdot 10^{-3}$ was not to rise above that of the regression start. The decision records the first part (messy failures 4 before, 1 after, the one left with another cause) and node date RMSE that moved both ways; it gives no separate RMSE comparison for that cell.
