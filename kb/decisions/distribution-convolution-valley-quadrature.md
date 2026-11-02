# Convolution samples below the FFT trust floor are computed by direct quadrature

**Type**: Intentional divergence from v0 (numerical method), approved by the user.

The time passes convolve two sampled time functions by FFT. An FFT result has an absolute round-off of about $10^{-16}$ of its peak, with either sign, so a sample far below the peak has no correct digits. v1 trusts the samples above the floor $10^{-13}$ of the peak (`CONV_TRUST_FRACTION`). When an operand has two modes, the result can fall below this floor between its modes and rise above it again: a **valley**. v1 computes every sample of the trusted range that lies at or below the floor by direct quadrature, the same trapezoid sum over the overlap of the operands that the FFT approximates ([distribution-convolution-trapezoid-rule.md](distribution-convolution-trapezoid-rule.md)). The tails outside the trusted range are computed by exponential tilting ([distribution-convolution-tails-by-exponential-tilting.md](distribution-convolution-tails-by-exponential-tilting.md)).

**v1 location**: `fn neg_log_with_quadrature_in_valleys()`, called from `fn reconstruct_neg_log_tails()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs); the quadrature of one output sample is `pub fn convolve_direct_at()` in [packages/treetime-ops/src/convolution.rs](../../packages/treetime-ops/src/convolution.rs).

**v0 location**: `NodeInterpolator.convolve_fft()` [packages/legacy/treetime/treetime/node_interpolator.py#L226-L265](../../packages/legacy/treetime/treetime/node_interpolator.py#L226-L265). v0 keeps only the samples above `fft_res.max() * 1e-13`, wherever they lie, and builds the result on that cropped grid with linear interpolation of the neg-log values. Across a valley, v0 therefore draws a straight line between the two trusted samples at the valley edges, both at about the trust floor: the valley density is the floor density, often many orders of magnitude above the true density.

## Why

Before this change, v1 took $-\ln$ of every sample between the first and the last trusted sample. Negative round-off samples in a valley gave NaN, and the time inference stopped ("Cannot re-peak a distribution ...: its values contain NaN"). Nextstrain norovirus builds (GII.4 genome, all genome, all RdRp, all VP1, and the example data all VP1) and `data/lassa/L/500` with the input root kept or with the clock filter off stopped this way.

Three treatments were measured: clamp each valley sample to the floor (the v0 behavior up to the slope between the two edge values), fill the valley in log space by continuing the slopes of both trusted edges, and direct quadrature. The valley holds less than $10^{-12}$ of the mass of the convolution in every measured case, but the valley values still change the results. In a later product, the valley of one message can lie under the peak of another factor, and the product then takes its shape from the valley values. Only direct quadrature gives correct valley values.

## Evidence

- **Closed form**: `test_convolution_of_bimodal_operand_with_deep_valley_matches_its_oracles` (`packages/treetime-distribution/src/distribution_ops/__tests__/test_edge_convolution.rs`) convolves an equal-weight mixture of two unit Gaussians at 20 and 80, truncated at $\pm 6$, with a unit Gaussian truncated at $\pm 6$, on a common spacing of 0.05. The oracles are the closed form of the truncated Gaussians (error functions) and an independent log-space trapezoid quadrature of the same samples:

  | quantity | bound | measured |
  | --- | --- | --- |
  | neg-log density at both modes, against the closed form | $10^{-10}$ | $3.3 \cdot 10^{-15}$ |
  | log total mass (trapezoid sum of the samples), against the closed form | $10^{-10}$ | $2.3 \cdot 10^{-11}$ |
  | samples within 6 of a mode, against the closed form | $10^{-6}$ | $8.7 \cdot 10^{-8}$ |
  | samples less than 20 nats below the peak, against the quadrature | $10^{-6}$ | $6.6 \cdot 10^{-8}$ |
  | valley samples more than 1 nat below the trust floor, against the quadrature | $10^{-9}$ | $3 \cdot 10^{-14}$ |

  At $t = 50$ the true valley lies 308.71 nats below the peak. Direct quadrature gives 308.71; the log-space fill gives 136.7, and the clamp gives 29.93, the trust floor. The closed form of the continuous functions differs from the quadrature of the samples by up to 0.065 nats in the deep valley, where the trapezoid rule meets the truncation of the kernel; this is the discretization error of the trapezoid rule
- **Production operands**: the 64 convolutions with valleys in the cases below (valleys 38 to 80 nats deep, operands of up to 97 000 and 30 000 samples) were captured and recomputed sample by sample by direct quadrature. The clamp overestimates the valley density by up to 50 nats (GII.4 genome), and the fill errs in both directions by up to 13.6 nats (all VP1). Direct quadrature matches the reference exactly, which here only checks the consistency of one formula; the closed-form test above is its independent evidence. No valley sample underflowed. The samples above the floor carry the FFT error of about $10^{-3}$ nats next to the floor, the same for every treatment
- **Runs**: release build, one thread, root dates and `log_lh_time` (the log marginal likelihood of node times) per treatment:

  | case | clamp | fill | direct |
  | --- | --- | --- | --- |
  | norovirus all genome | 1930.0, -345740.4 | 1948.2, -345679.7 | 1951.9, -345654.1 |
  | norovirus all RdRp | 1535.8, -213820.5 | 1525.9, -213820.9 | 1538.8, -213819.8 |
  | norovirus all VP1 | 1722.3, -354017.2 | 1725.4, -354011.4 | 1722.1, -354010.7 |
  | lassa L 500, Nextstrain pipeline config | 1551.10, -6782.9 | 1551.24, -6782.2 | 1551.20, -6781.1 |

  The norovirus GII.4 genome and example data builds and three `lassa/L/500` settings (input root kept with the coalescent prior, input root or midpoint root with the clock filter off) give the same dates under all three treatments. With the midpoint root, the coalescent prior and an estimated clock rate, the rate after two iterations is $6.91 \cdot 10^{-4}$ with the clamp and $7.90 \cdot 10^{-4}$ with direct quadrature; that run has not converged after its two iterations, so its spread is not a clean comparison
- **Cost**: direct quadrature costs one overlap sum per valley sample. The run time grows by 0 to 30% (largest: 10.3 s to 13.4 s on `lassa/L/500` with the input root kept and the coalescent prior), and the peak memory does not change
- **v0** dates the same builds 11 years (lassa pipeline config) to about 400 years (norovirus all VP1) away from every treatment, because of other differences between v0 and v1, so v0 cannot select a treatment

## Limits

- A valley sample whose quadrature underflows (more than about 708 nats below the peak) is set to the smallest normal double; no measured case reaches this depth
- The trusted range ends at the outermost samples above the floor. Below-floor samples outside it belong to the tails, which exponential tilting computes

## Related

- [distribution-convolution-trapezoid-rule.md](distribution-convolution-trapezoid-rule.md): the quadrature that both the FFT and the direct sum compute
- [kb/issues/N-distribution-mixed-nan-policy-undecided.md](../issues/N-distribution-mixed-nan-policy-undecided.md): the policy for NaN in sampled values
