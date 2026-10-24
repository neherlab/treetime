# The timetree start rate falls back to the Gaussian rate below a divergence floor

With an estimated clock rate, `timetree` starts from the slope of the root-to-tip regression. It falls back to the
substitution rate of the Gaussian timetree when that slope is at most

$$
\mu_\text{floor} = 0.1\,\frac{\bar d}{t_\text{max} - t_\text{min}}
$$

where $\bar d$ is the mean root-to-tip divergence of the dated tips that are not clock-filter outliers, for the root
the fit uses, and $t_\text{max} - t_\text{min}$ their date span. A slope at or below it would have the divergence
accumulate over at least ten sampling windows. The floor is 0 when the dates do not vary, which keeps the earlier
rule of falling back only for a slope at or below zero.

**Type**: Scientific default, new in v1. v0 starts from the regression slope.

**v1 location**: `fn ClockRerootResult::signal_rate_floor()`
[packages/treetime/src/clock/clock_regression.rs#L176](../../packages/treetime/src/clock/clock_regression.rs#L176)
(fraction `SIGNAL_RATE_FLOOR_FRACTION` [#L28](../../packages/treetime/src/clock/clock_regression.rs#L28)),
`fn ClockRerootResult::into_clock_fit_or()`
[#L157](../../packages/treetime/src/clock/clock_regression.rs#L157), called with the floor by
`fn estimate_initial_clock()`
[packages/treetime/src/timetree/pipeline.rs#L278-L305](../../packages/treetime/src/timetree/pipeline.rs#L278-L305);
the fallback rate is `fn gaussian_start_rate()`
[packages/treetime/src/timetree/init.rs#L186-L210](../../packages/treetime/src/timetree/init.rs#L186-L210).

## Why

On a tree without temporal signal the regression slope scatters around zero. A slightly negative slope fell back to
the Gaussian rate, a slightly positive one was kept, and no round recovers a start orders of magnitude too low: the
first round moves the rate at most 2-fold ([timetree-clock-rate-by-marginal-likelihood.md](timetree-clock-rate-by-marginal-likelihood.md#rate-window)).
On messy validation data (`validation/degrade.py`; `data_messy/const/rep03`, $L=10^3$, true rate $3\cdot10^{-3}$,
`--keep-root --clock-filter 0`) the slope was $2.29\cdot10^{-5}$: without a coalescent prior the run stopped at the
grid limit, and with `--coalescent opt` or `skyline` it panicked in the date output after dating nodes before year 0
([H-timetree-negative-numdate-panics-in-date-string-conversion.md](../issues/H-timetree-negative-numdate-panics-in-date-string-conversion.md)). The
floor there is $2.57\cdot10^{-3}$, the Gaussian rate $3.95\cdot10^{-3}$, and the runs end at 0.96 (no prior), 0.98
(`opt`) and 1.01 (`skyline`) of the true rate.

The start is only a start: the refinement loop optimizes the rate further, so the floor needs to catch a hopeless
start, not to be exact.

The floor applies only to the start. Inside the refinement loop the regression runs on clock lengths (time times the
current rate), so the floor would follow the current rate; the in-loop refits keep the threshold of zero.

## Validation

Timetree `core` set, estimated rate, the six variants `{keep,reroot}-{none,opt,skyline}-est`, 400 tips (240 sets
each), `origin/rust` against this change:

| data | runs | fallback fired | failed before / after |
| --- | --- | --- | --- |
| clean (`data`) | 1440 | 14 | 0 / 0 |
| messy (`data_messy`) | 1440 | 15 | 4 / 1 |

The floor fired only in `keep-*` runs; with `--reroot` the root search finds a root with a positive slope. Where it
fired with a coalescent prior, the final rate moved toward the truth in 17 of 18 runs, for example on `const/rep02`,
$L=10^3$, $\mu=3\cdot10^{-4}$ ($\ln$ rate error $-0.61 \to +0.45$ with `opt`, $-0.58 \to +0.45$ with `skyline`) and
on `const/rep05`, $L=10^4$, $\mu=3\cdot10^{-4}$ ($-0.28 \to +0.10$); the exception is `keep-skyline-est` on messy `const/rep02`, $L=10^3$, $\mu=3\cdot10^{-4}$
($-0.36 \to +0.44$). The node date RMSE moved both ways. The messy
failure left is `keep-none-est` on `const/rep03`, $L=10^3$, $\mu=10^{-3}$: its start is the Gaussian rate before and
after (negative slope), and the in-loop regression of undated tips then runs the rate down
([M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md](../issues/M-timetree-clock-update-regresses-undated-tips-on-inferred-dates.md)).

## Limitation: deep trees with a short sampling window

$\bar d / (t_\text{max} - t_\text{min}) \approx \mu \cdot \text{depth}/\text{span}$, so the floor lies below the true
rate only while the root is less than about ten sampling spans deep. On the validation trees it is 0.2 (median) to 0.86
of the true rate. Keeping only the last 0.3 years of samples of 20 trees of `data_n30` (120 tips, root 18 to 100
spans deep), the floor lay at 1.8 to 10 times the true rate and fired on every set, although the regression slope was
$-0.7$ to $2.6$ times the true rate. The Gaussian substitution rate is biased high there (1.4 to 2.6 times the true
rate): its conversion $\mu = \eta D/(D+M+N)$ assumes one hazard for the branching and sampling events of all branches,
and a dense recent sample has a much higher sampling hazard on terminal branches than the deep internal branches have
(Section 15.4 of the maintainers' derivation of the Gaussian timetree). Without a coalescent prior the final rates were
unchanged; with `opt` the median $\lvert\ln\rvert$ final rate error rose from 0.115 to 0.285, with no failures. A
two-class hazard in that conversion is the documented remedy and is not implemented.
