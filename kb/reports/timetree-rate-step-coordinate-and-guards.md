# Coordinate and guards of the marginal clock rate step

Date: 2026-10-23

With a coalescent prior and an estimated rate, each refinement round moves the clock rate toward the maximum of the log marginal likelihood of the node times ($\log Z$, `log_lh_time`), with a parabola through $\log Z$ at the rate and at the rate divided and multiplied by 1.25 ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)). The parabola is fitted in log rate, and three guards bound the step: the target is capped at a factor 2, from the second round only a maximum between the probes counts (`RateWindow::InteriorMaximum`), and the rate moves 2/3 of the way to the target. The guards were chosen while $\log Z$ carried a first-order convolution error that tilted it toward slow rates ([kb/reports/timetree-low-signal-rate-bias-from-convolution-edges.md](timetree-low-signal-rate-bias-from-convolution-edges.md)). This report asks whether the guards are still needed now that the trapezoid rule removed that error, and whether a parabola in another coordinate of the rate fits $\log Z$ better.

## Summary

- In log rate, the fixed point of the step is not the maximum of $\log Z$: it lies $h^2 f'''/(6 f'')$ below it ($h = \ln 1.25$), which for the Poisson shape of $\log Z$ is $-h^2/6 = -0.8\%$ of the rate. The validation runs show this bias at high signal ($-1.1\%$ to $-1.3\%$, standard error $0.1\%$ to $0.2\%$), and a parabola in the cube root or square root of the rate removes most of it
- On profiles of $\log Z$, the cube root is the most nearly quadratic coordinate near the maximum, the square root next; log rate is skewed the other way, and the rate itself most. This is the order that a Poisson likelihood predicts
- Dropping the interior window (following the parabola maximum in every round, within the 2-fold cap) changes no run at high signal and lowers the median rate error at low signal ($L=10^3$, $\mu=3\cdot10^{-4}$) from 0.30 to 0.25
- A cube-root parabola with the full step gives the lowest rate error overall (median $\lvert\ln\rvert$ error 0.026 against 0.033 on 400 tips, 0.022 against 0.029 on 1200 tips). At the lowest signal on 400 tips more runs end beyond 2-fold of the true rate (15% against 11%); on 1200 tips fewer do (12.5% against 28%)
- At low signal $\log Z$ jumps by several nats between neighboring rates, which limits every three-probe step whatever its coordinate
- No variant of the step recovers a start rate that is far off: on the case of [kb/issues/M-timetree-start-rate-accepts-near-zero-regression-slope.md](../issues/M-timetree-start-rate-accepts-near-zero-regression-slope.md) the coalescent timescale fits the stretched tree, and $\log Z$ under that prior peaks near the wrong rate

## Fixed point of the three-probe step

With probes at $x-h$, $x$ and $x+h$ in log rate $x$, the parabola maximum lies at $x + h\,(f_- - f_+)/(2(f_+ + f_- - 2f_0))$. The step leaves the rate in place when $f_+ = f_-$. From $f_+ - f_- = 2hf'(x) + h^3 f'''(x)/3 + O(h^5)$, that happens at $f'(x) = -h^2 f'''/6$, that is at

$$x - x^* = -\frac{h^2 f'''}{6 f''}$$

from the maximum $x^*$. A partial step converges geometrically to this point, not to $x^*$. If $\log Z$ were quadratic in some coordinate $u$ of the rate, the parabola in $u$ through the three probes would be exact, and its fixed point the maximum.

For one Poisson term $k \ln\lambda - \lambda$ in the coordinate $\phi = \lambda^{1/a}$, the third derivative at the maximum, standardized by the curvature, is $(2 - (a-1)(a-2))/(a^2\sqrt k)$: $2/\sqrt k$ for the rate, $0.5/\sqrt k$ for the square root, 0 for the cube root, and $-1/\sqrt k$ for log rate ($a \to 0$, where $f'''/f'' = 1$ and the offset is $-h^2/6$). With node times fixed, the substitution part of $\log Z$ is a sum of such terms in the rate; the integration over node times and the coalescent prior add to it.

## Profiles of log Z

Profiles from the ignored test `research_rate_profile_shape` (branch `research/rate-step`): validation layout (10 tips at each of $\tau = 0, 0.1, \ldots, 3.9$, constant $T_c = 10$, random GTR), 5 seeds per $(L, \mu)$, $\log Z$ at $x = \ln(\text{rate}/\mu) \in [-1.6, 1.6]$ in steps of 0.1, under the coalescent prior frozen at the true node times with the branch likelihoods at the true rate (`FrozenRateState`). Analysis: `test_scripts/timetree_rate_step_shape.py` (branch `research/rate-step`). Medians over seeds:

| $L$ | $\mu$ | $x^*$ | $f''$ in $x$ | std. $f'''$: rate | sqrt | cbrt | log |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $10^3$ | $3\cdot10^{-4}$ | 0.217 | -960 | -0.04 | -0.11 | -0.13 | -0.17 |
| $10^3$ | $10^{-3}$ | 0.030 | -486 | 0.99 | 0.93 | 0.90 | 0.86 |
| $10^3$ | $3\cdot10^{-3}$ | 0.005 | -507 | 0.069 | -0.003 | -0.027 | -0.076 |
| $10^3$ | $10^{-2}$ | -0.006 | -1185 | 0.055 | 0.010 | -0.005 | -0.035 |
| $10^4$ | $3\cdot10^{-4}$ | 0.003 | -482 | 0.209 | 0.138 | 0.114 | 0.067 |
| $10^4$ | $10^{-3}$ | -0.001 | -1171 | 0.070 | 0.025 | 0.011 | -0.019 |
| $10^4$ | $3\cdot10^{-3}$ | 0.001 | -2989 | 0.034 | 0.008 | -0.001 | -0.019 |
| $10^4$ | $10^{-2}$ | -0.007 | -9222 | 0.018 | 0.002 | -0.003 | -0.013 |

The maximum lies at the true rate within 1% except at the lowest signal. Where the profile is smooth, the standardized third derivative follows the predicted order rate > sqrt > cbrt $\approx$ 0 > log.

Error $\lvert x_1 - x^*\rvert$ after one parabola step started 2-fold away ($x_0 = x^* \pm 0.7$), median over seeds:

| $L$ | $\mu$ | start | rate | sqrt | cbrt | log |
| --- | --- | --- | --- | --- | --- | --- |
| $10^3$ | $3\cdot10^{-3}$ | below | 0.29 | 0.14 | 0.06 | 0.22 |
| | | above | no maximum | 0.23 | 0.05 | 0.15 |
| $10^4$ | $10^{-3}$ | below | 0.29 | 0.13 | 0.05 | 0.25 |
| | | above | 1.31 | 0.02 | 0.09 | 0.23 |
| $10^4$ | $3\cdot10^{-3}$ | below | 0.26 | 0.07 | 0.04 | 0.46 |
| | | above | 0.05 | 0.23 | 0.28 | 0.34 |

Repeated full steps without guards converge (step below 0.01) in the cube root from up to 4-fold away at every $(L,\mu)$ but the lowest signal, in a median of 3 to 5 steps (7 from 4-fold above at $L=10^3$, $\mu=10^{-3}$); in log rate, starts 4-fold below overshoot out of the profile in some sets, and in the rate, starts above the maximum often find no maximum. Over wider windows ($\lvert x - x^*\rvert \le 0.7$ and 1.2) the cube root fits best except at $L=10^4$, $\mu \ge 3\cdot10^{-3}$, where the curvature in $x$ grows faster than $e^x$ above the maximum and the square root fits best.

At $L=10^3$, $\mu=3\cdot10^{-4}$, second differences of $\log Z$ at spacing 0.1 scatter by $\pm 5$ to $\pm 25$ nats around a smooth trend, against a smooth progression at $L=10^4$. Between neighboring rates $\log Z$ jumps by several nats (for example $-19.8$ to $-3.1$ between $x = -0.4$ and $-0.2$ in one seed). A rescaled rate moves the time grids, so this is likely the grid sensitivity of [kb/issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md](../issues/M-timetree-node-dates-move-when-a-time-grid-gains-one-point.md). No coordinate helps there.

## Validation runs

Timetree `core` set, 10 tips per sample time, clean data, the four variants with a coalescent prior and an estimated rate (`{keep,reroot}-{opt,skyline}-est`), 240 sets each, 960 runs per row; research build of branch `research/rate-step` at `0c2fe5ab5`. The base row reproduces the runs of the unmodified build exactly. "No window": `RateWindow::FollowSlope` in every round, with the 2-fold cap. Analysis: `test_scripts/timetree_rate_step_compare.py` (branch `research/rate-step`).

Median $\lvert\ln(\text{rate}/\mu)\rvert$:

| step | $10^3$, $3\cdot10^{-4}$ | $10^3$, $10^{-3}$ | $10^3$, $3\cdot10^{-3}$ | $10^3$, $10^{-2}$ | $10^4$, $3\cdot10^{-4}$ | $10^4$, $10^{-3}$ | $10^4$, $3\cdot10^{-3}$ | $10^4$, $10^{-2}$ | all |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| log, window, 2/3 (current) | 0.302 | 0.088 | 0.041 | 0.032 | 0.049 | 0.024 | 0.016 | 0.012 | 0.033 |
| log, no window, 2/3 | 0.250 | 0.075 | 0.041 | 0.032 | 0.048 | 0.024 | 0.016 | 0.012 | 0.032 |
| sqrt, no window, 2/3 | 0.253 | 0.077 | 0.034 | 0.029 | 0.049 | 0.024 | 0.012 | 0.007 | 0.029 |
| cbrt, no window, 2/3 | 0.245 | 0.076 | 0.036 | 0.031 | 0.048 | 0.023 | 0.013 | 0.008 | 0.030 |
| cbrt, no window, full step | 0.275 | 0.076 | 0.031 | 0.030 | 0.034 | 0.020 | 0.014 | 0.007 | 0.026 |

Mean signed $\ln(\text{rate}/\mu)$ at high signal (standard error 0.001 to 0.002):

| step | $10^4$, $10^{-3}$ | $10^4$, $3\cdot10^{-3}$ | $10^4$, $10^{-2}$ |
| --- | --- | --- | --- |
| log, window, 2/3 (current) | -0.0051 | -0.0112 | -0.0129 |
| log, no window, 2/3 | -0.0051 | -0.0112 | -0.0129 |
| sqrt, no window, 2/3 | +0.0021 | -0.0014 | -0.0032 |
| cbrt, no window, 2/3 | -0.0003 | -0.0046 | -0.0064 |
| cbrt, no window, full step | -0.0013 | -0.0038 | -0.0058 |

The cube root removes about 0.0065 of the bias, close to the predicted $h^2/6 = 0.0083$; the square root removes more. The bias left in the cube root has another cause.

Other measures, all 960 runs per row:

| step | sets with smaller / larger rate error than current | within 10% of $\mu$ at $10^3$, $3\cdot10^{-3}$ | beyond 2-fold at $10^3$, $3\cdot10^{-4}$ | median date RMSE (years) | median root error (years) |
| --- | --- | --- | --- | --- | --- |
| log, window, 2/3 (current) | | 0.81 | 0.108 | 0.118 | 0.569 |
| log, no window, 2/3 | 0.06 / 0.03 | 0.81 | 0.058 | 0.118 | 0.569 |
| sqrt, no window, 2/3 | 0.61 / 0.39 | 0.82 | 0.067 | 0.118 | 0.564 |
| cbrt, no window, 2/3 | 0.61 / 0.39 | 0.80 | 0.058 | 0.118 | 0.570 |
| cbrt, no window, full step | 0.69 / 0.31 | 0.93 | 0.150 | 0.117 | 0.533 |

No run failed in any variant. The node date RMSE barely moves: at these rate errors the dates are limited by other errors.

The same comparison with 30 tips per sample time (`data_n30`, 1200 tips; no run failed):

| step | median $\lvert\ln\rvert$ rate error, all | signed bias at $10^4$, $3\cdot10^{-3}$ / $10^{-2}$ | beyond 2-fold at $10^3$, $3\cdot10^{-4}$ | smaller / larger rate error than current | median date RMSE (years) |
| --- | --- | --- | --- | --- | --- |
| log, window, 2/3 (current) | 0.029 | -0.0130 / -0.0111 | 0.283 | | 0.076 |
| sqrt, no window, 2/3 | 0.029 | -0.0032 / -0.0009 | 0.150 | 0.57 / 0.43 | 0.077 |
| cbrt, no window, 2/3 | 0.028 | -0.0064 / -0.0043 | 0.150 | 0.58 / 0.42 | 0.076 |
| cbrt, no window, full step | 0.022 | -0.0041 / -0.0033 | 0.125 | 0.63 / 0.37 | 0.075 |

The high-signal bias and its removal repeat. At the lowest signal, dropping the window halves the share of runs beyond 2-fold, and here the full step does not raise it.

## Far-off start

On messy data, `keep-opt-est` on `const/rep03` at $L=10^3$, $\mu=3\cdot10^{-3}$ starts at $2.29\cdot10^{-5}$ ([kb/issues/M-timetree-start-rate-accepts-near-zero-regression-slope.md](../issues/M-timetree-start-rate-accepts-near-zero-regression-slope.md)). Every variant above, also with the cap raised to 10-fold and the full step, ends between $4.1\cdot10^{-5}$ and $4.6\cdot10^{-5}$ and panics in the date output. With the cube root, full step and 10-fold cap, the two rounds move the rate by factors 1.30 and 1.53, each to the parabola maximum: the coalescent timescale fitted to the stretched tree is $T_c \approx 5500$ years (true: 10), and $\log Z$ under that prior peaks near the current rate. The loop sits at a self-consistent wrong state, which no step rule leaves; the start rate has to be right.

## Options

- **Coordinate:** fit the parabola in the cube root or the square root of the rate. Both remove most of the high-signal bias; the cube root fits better near the maximum, the square root over wide windows at the highest signal and removes more of the bias in the runs
- **Window:** drop `RateWindow::InteriorMaximum`. It changes nothing where the signal is clear and was the cause of a worse rate at low signal
- **Step fraction:** the full step converges in fewer rounds at medium and high signal; at the lowest signal it ends more runs beyond 2-fold of the true rate than 2/3 on 400 tips, and fewer on 1200 tips
- **Low-signal noise:** a least-squares parabola over more probes, or wider probes, would average the grid noise of $\log Z$; this is the option already listed for the flat top of $\log Z$ in [kb/reports/timetree-low-signal-rate-bias-from-convolution-edges.md](timetree-low-signal-rate-bias-from-convolution-edges.md)

## Reproduction

Branch `research/rate-step` (local), commits `0c2fe5ab5` (research build and profile test) and `87c96797b` (scripts).

- Profiles: `cargo test --release -p treetime --lib -- --ignored --nocapture --exact test_utils::__tests__::test_rate_profile::tests::research_rate_profile_shape > profile.tsv` (18 minutes on 4 threads), then keep the TSV lines and run `test_scripts/timetree_rate_step_shape.py`
- Validation: environment of the research build, per variant: unset (current), `TREETIME_RATE_GUARD=none` (no window), with `TREETIME_RATE_COORD=sqrt` or `cbrt`, and `TREETIME_RATE_STEP=1` (full step); commands in the docstring of `test_scripts/timetree_rate_step_compare.py`

## Follow-up

The step now keeps the coordinate (log rate) and the gate, and halves the probe half-width after each step to a maximum between the probes, down to 1/8 of the start. That removes most of the high-signal offset without the parabola in the cube root, and without the gate the rate drifts at low signal over six rounds: [kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#shrinking-probe-window).
