# Convolution valley quadrature turns an exact zero into a finite density

The function-function convolution computes every sample of its trusted range that lies at or below the FFT trust floor by direct quadrature ([kb/decisions/distribution-convolution-valley-quadrature.md](../decisions/distribution-convolution-valley-quadrature.md)). `fn neg_log_with_quadrature_in_valleys()` clamps the quadrature result to the smallest normal double before it takes the negative logarithm: `convolve_direct_at(dx, pa, pb, i)?.max(f64::MIN_POSITIVE)` ([packages/treetime-distribution/src/distribution_ops/convolve.rs#L425-L444](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L425-L444), clamp at [#L439](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L439)). A quadrature result of exactly zero therefore becomes a finite cost of $-\ln(2.2 \cdot 10^{-308}) \approx 708.4$ nats above the peak-normalized operands, where the correct value is $+\infty$.

## When the quadrature returns zero

The quadrature `pub fn convolve_direct_at()` ([packages/treetime-ops/src/convolution.rs#L88-L108](../../packages/treetime-ops/src/convolution.rs#L88-L108)) sums the products of the plain operand values `pa` and `pb`. The sum is exactly zero when, for every pair of samples that overlap at the output sample, one of the two values is zero. Two different causes give a zero plain value, and the quadrature cannot tell them apart:

- **Support gap**: a `DistributionFunction` can hold `+inf` neg-log values inside its grid, so one operand can have two finite-support components with an impossible region between them (`fn narrowest_support_component()` and `fn convolution_range_function()` in the same file handle such interior infinite values). If the gap is wider than the support of the other operand, the true convolution is zero between the two result components. The trusted range runs from the first to the last sample above the floor ([#L371-L376](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L371-L376)), so the gap lies inside it and every gap sample gets about 708 nats instead of $+\infty$. A later product with a factor whose mass lies in the gap then puts probability into a region that the convolution excludes
- **Underflow**: `PeakNormalizedOperand::new()` computes the plain values as $e^{-r}$ of the neg-log values relative to the operand peak ([#L519](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L519)). Samples more than about 745 nats below the peak underflow to zero. An operand with two modes and a valley deeper than about 745 nats has a run of zero plain values, and the convolution valley can then be zero in floating point while the true density is positive. The clamp then gives 708 nats where the true cost is larger than 745 nats, so the valley also gets too much density

The decision file lists the underflow case under "Limits" and notes that no measured case reaches this depth (the measured production valleys lie 38 to 80 nats below the peak). It does not cover the support gap, where the correct result is an exact zero and not a small positive number.

> [!IMPORTANT]
> **Investigation required.** It is not shown that either cause is reachable from `timetree` inputs. The time messages of the backward and forward passes have one support interval in exact arithmetic: leaf date constraints are points or ranges, the branch-length likelihoods (`packages/treetime/src/timetree/inference/input_branch_length.rs`) are infinite only at the end of their grid, products intersect the supports, and convolutions of intervals are intervals; the support of the cavity division of the forward pass is not checked. No code path that writes an interior `+inf` into a time message was found. Collect: (a) a time message of a real or simulated run with an interior `+inf` neg-log value, or a proof that none exists, and (b) the largest valley depth of a convolution operand on the example datasets and the Nextstrain cases of [kb/decisions/distribution-convolution-valley-quadrature.md](../decisions/distribution-convolution-valley-quadrature.md), against the underflow depth of about 745 nats.

## Fix direction

- Keep an exact zero of the quadrature as $+\infty$ in neg-log space, so a support gap stays impossible
- Compute the valley quadrature in log space (log-sum-exp over the products of the overlap, from the relative neg-log values instead of the plain values), so that underflow can be told apart from an exact zero and a deep valley keeps its true depth
- Add a convolution test with an exact interior support gap: an operand with two finite-support components separated by `+inf` values, convolved with an operand whose support is narrower than the gap; the result must be `+inf` at every sample of the gap. Its oracle is the support arithmetic of the two operands
- Update the "Limits" section of the decision file when the fix changes the underflow behavior

## Related issues

- [N-distribution-plain-division-fixed-floor.md](N-distribution-plain-division-fixed-floor.md): the Plain division clamps zero divisors to a fixed floor, the same pattern of a clamp that turns zero into a finite value
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md): the policy for NaN in sampled values
- [M-distribution-convolution-drops-operand-tails.md](M-distribution-convolution-drops-operand-tails.md): the tails of the convolution operands
