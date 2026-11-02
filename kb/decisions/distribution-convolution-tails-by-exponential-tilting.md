# Convolution tails below the FFT trust floor are computed by exponential tilting

**Type**: Intentional divergence from v0 (numerical method), approved by the user.

The time passes convolve two sampled time functions by FFT. An FFT result has an absolute round-off of about $10^{-16}$ of its peak, so v1 trusts only the samples above the floor $10^{-13}$ of the peak (`CONV_TRUST_FRACTION`), about 30 nats below it. The samples beyond the outermost trusted samples are the **tails** of the convolution. v1 computes them by exponential tilting: it multiplies both operands by $e^{\theta t}$, convolves them again by FFT, and divides the result by $e^{\theta t}$. With $\theta$ equal to the slope of the resolved neg-log values at the end of the tail, the tilted convolution is flat there, so the next samples rise above the trust floor of the tilted result. v1 repeats this with the slope at the new end until the tail reaches the end of the convolution or a tilt adds no sample. At an end where both operands end hard, the result ends at the outermost resolved sample, so the samples that no tilt reaches get zero density. At an end where an operand continues with a soft tail law, they get a straight-line extension of the outermost resolved samples. The valley samples inside the trusted range are computed by direct quadrature ([distribution-convolution-valley-quadrature.md](distribution-convolution-valley-quadrature.md)).

**v1 location**: `fn neg_log_at_ends_by_tilting()` and `struct TiltedConvolution`, called from `fn reconstruct_neg_log_tails()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs).

**v0 location**: `NodeInterpolator.convolve_fft()` [packages/legacy/treetime/treetime/node_interpolator.py#L226-L265](../../packages/legacy/treetime/treetime/node_interpolator.py#L226-L265). v0 keeps the samples above `fft_res.max() * 1e-13` and extends the neg-log values beyond them by straight lines through the outermost trusted sample, with the slope over the last 3 trusted samples. v1 used the same straight lines before this change.

## Method

Let $f$ and $g$ be the operand densities on a common spacing $\Delta t$, with samples $f_i$ and $g_j$. The trapezoid sum of the convolution ([distribution-convolution-trapezoid-rule.md](distribution-convolution-trapezoid-rule.md)) at output sample $n$ is

$$c_n = \Delta t \sum_{i} w_{n,i}\, f_i\, g_{n-i}$$

where:

- $c_n$ -- convolution sample $n$
- $w_{n,i}$ -- trapezoid weight, which depends only on the indices $n$ and $i$
- $f_i$, $g_j$ -- operand samples at $t_i = i \Delta t$ and $t_j = j \Delta t$ from the start of each grid

The tilted operands $\tilde f_i = f_i\, e^{\theta t_i}$ and $\tilde g_j = g_j\, e^{\theta t_j}$ give $\tilde c_n = c_n\, e^{\theta t_n}$ exactly, because $t_i + t_{n-i} = t_n$ for every term. In neg-log values, $-\ln c_n = -\ln \tilde c_n + \theta t_n$. The tilt changes which samples lie above the trust floor and not the quadrature. v1 computes the tilted operands from the neg-log values, so operand samples that underflow in the untilted operands keep their values.

## Why

A chain of internal branches without substitutions is a chain of zero-length branches. Along such a chain, the subtree message of a node lies deep in the young tail of the convolution that sends it to the parent: for node 100 of a 200-tip caterpillar, 60.5 nats below the peak. The straight-line extension is far too heavy there, the tail then holds the mass of the message, and the posteriors come out too young and too narrow. The error does not depend on the grid: refining the branch grid leaves it in place.

## Evidence

- **Closed form**: a caterpillar of $n$ tips that share one date, with no substitution on any branch, has $\log Z = -\sum_{m=2}^{n} \ln(\lambda m)$, where $Z$ is the marginal likelihood of the node times and $\lambda$ the expected substitutions per year on one branch ([kb/issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](../issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). Error of $\log Z$ in nats on the 200-tip caterpillar ($\lambda = 1.0556$ per year), against the number of points of the branch grid:

  | branch grid points | straight line (v0 rule) | tilting |
  | ---: | ---: | ---: |
  | 300 (default) | 171.5 | 73.9 |
  | 599 | 150.2 | 27.7 |
  | 1197 | 145.0 | 8.13 |
  | 2393 | 143.2 | 2.14 |
  | 4785 | 142.8 | 0.546 |

  With straight lines the error stays near 143 nats. With tilting it falls by 2.7, 3.4, 3.79 and 3.93 times per halving of the spacing and approaches the second order of the trapezoid rule. On 20 and 50 tips, the messages stay above the trust floor where the posteriors need them, and both rules give the same errors (0.116 and 1.89 nats at 300 points), which fall by 3.75 to 3.96 times per halving. `test_zero_length_ladder_log_marginal_likelihood_converges_at_second_order_to_the_exact_value` asserts a fall of at least 3 times per halving at three spacings. The 400-tip caterpillar stopped at the grid point limit with straight lines and dates every node with tilting
- **Fine grid**: on the smoke cases of `timetree`, release builds at 300 and 1200 grid points (`GRID_POINTS`), maximum and mean distance of the internal node dates from the tilting run at 1200 points, in years, and `log_lh_time`:

  | case | straight line, 300 | tilting, 300 | straight line, 1200 | `log_lh_time`: straight line, tilting at 300; tilting at 1200 |
  | --- | --- | --- | --- | --- |
  | `lassa/L/50/basic` | 10.2, 0.89 | 0.084, 0.025 | 10.1, 0.88 | -272.94, -287.18; -287.18 |
  | `ebola/20/coal-opt` | 0.18, 0.022 | 0.0009, 0.0005 | 0.0002, 0.0000 | -72.02, -75.14; -75.13 |
  | `lassa/L/20/reroot-midpoint` | 1.93, 0.34 | 0.85, 0.14 | 1.63, 0.18 | -91.894, -91.897; -91.888 |
  | `rsv/a/100/basic` | 0.090, 0.013 | 0.066, 0.008 | 0.0016, 0.0003 | -160.607, -160.601; -159.877 |
  | `dengue/100/basic` | 0.044, 0.006 | 0.042, 0.006 | 0.0066, 0.0019 | -189.359, -189.359; -189.361 |
  | `tb/20/basic` | 0.51, 0.069 | 0.38, 0.049 | 0.0000, 0.0000 | 24.783, 24.783; 24.797 |
  | `dengue/20/basic` | 0.22, 0.058 | 0.48, 0.12 | 0.0037, 0.0010 | -228.055, -228.073; -228.417 |
  | `zika/20/clock-rate-1e-6` | 2.30, 0.76 | 1.61, 0.85 | 0.71, 0.22 | 68.073, 68.072; 68.081 |

  On `lassa/L/50/basic` the straight lines move dates by up to 10 years and the log marginal likelihood by 14 nats at every grid size, while tilting agrees with its fine-grid run. Where both rules agree at 1200 points (`dengue/20/basic`, `tb/20/basic`), the tails do not matter, and both differ from the fine grid by the discretization error of the default grid. Of the 52 `timetree` smoke cases that change, the others move dates by at most 0.033 years, except `tb/100/basic` and two `ebola/362` cases, whose clock filter and polytomy resolution take other paths; there, both rules differ from the fine grid by more than 1 year
- **Exact posteriors of a four-tip tree** (`test_marginal_exact`): the CDF differences to the exact posterior without a coalescent prior change by $-4 \cdot 10^{-5}$ to $+4 \cdot 10^{-6}$, and those of the mean-field coalescent to the exact Kingman posterior by $-3 \cdot 10^{-6}$ to $+4 \cdot 10^{-5}$, against differences of $2 \cdot 10^{-5}$ to $1.6 \cdot 10^{-2}$; the mean-field differences hold the model error of the mean-field coalescent
- **Production runs**: Nextstrain builds whose input trees hold chains of up to 288 (rubella E1), 421 (oropouche S), 428 (H5N1 cattle-outbreak HA), and 511 (PB2) zero-length internal branches. With the straight lines and without the pre-loop short-branch collapse, all four stopped at the grid point limit (1 003 324 to 1 026 440 points requested). With tilting, release build, fixed clock rates of the workflows:

  | build | without the collapse | with the collapse (default) |
  | --- | ---: | ---: |
  | rubella E1 | 98 s | 28 s |
  | oropouche S | 11 s | 3.1 s |
  | HA | 58 s | 21 s |
  | PB2 | 122 s | 24 s |

- **Cost**: a convolution whose tails lie above the trust floor runs no tilt. Each tilt costs one FFT convolution of the full operands. On `data/mpox/clade-ii/1000` the run time grows from 22.2 to 23.0 seconds and on `data/flu/h3n2/500` from 0.37 to 0.43 seconds (release build, two runs each); the peak memory does not change

## Limits

- Each tilt reads the samples of its result up to its own trust floor, so the samples next to that floor carry the FFT round-off of about $10^{-3}$ nats, as the untilted samples next to the untilted floor do
- The tilting removes the tail error. The second-order discretization error remains, and it grows fast along chains: 73.9 nats on the 200-tip caterpillar at the default grid ([kb/issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](../issues/H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md))

## Related

- [distribution-convolution-valley-quadrature.md](distribution-convolution-valley-quadrature.md): the samples below the trust floor inside the trusted range
- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md): the tail laws of the messages and the error sources of the log marginal likelihood
