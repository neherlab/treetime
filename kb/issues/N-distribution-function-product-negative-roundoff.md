# Function-function product yields tiny negative amplitudes from interpolation rounding

Multiplying two non-negative `Distribution::Function` operands can produce grid amplitudes that are slightly negative (order `-1e-26`), even though the exact product of two non-negative functions is non-negative everywhere. The negatives are pure floating-point rounding in the piecewise-linear interpolation, at grid points near the support tail where amplitudes are already tiny (order `1e-13`).

## Mechanism

`multiply_function_function` [packages/treetime-distribution/src/distribution_ops/multiply.rs#L150](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L150) resamples onto a fresh uniform grid over the support intersection and evaluates `a.interp(t) * b.interp(t)` at each point. `GridFn::interpolate_at` [packages/treetime-grid/src/grid_fn.rs#L338](../../packages/treetime-grid/src/grid_fn.rs#L338) computes `y0 + t * (y1 - y0)`. For a decreasing tail (`y1 < y0`, both near `1e-13`) with `t` at the interval edge, catastrophic cancellation and `t = (q - x0) / dx` rounding can push the interpolated value a few ulps below zero. The product of that tiny negative with the other operand's positive value yields a `~-1e-26` amplitude, which `normalize()` preserves (it only rescales by the maximum).

## Evidence

Traced on `timetree --tree data/ebola/20/tree.nwk --aln data/ebola/20/aln.fasta.xz --metadata data/ebola/20/metadata.tsv --coalescent 1.0`, instrumenting the backward pass to report the minimum amplitude at each stage:

- `stage=parent_message` never reports a negative: the child messages entering the fold are non-negative.
- `stage=product` reports `min = -4.79e-26 ... -8.69e-26` at the same internal node while folding one child: the negative is created by `distribution_multiplication(...).normalize()`.
- `stage=convolution` shows a negative `node_time_min` feeding in but a positive `conv_min` out: the convolution neither creates nor propagates the negative.

## Impact

Benign for current consumers. The magnitude is 25+ orders below the unit peak. The one operation that previously mishandled it, `distribution_apply_neg_log_weight`, took `ln(amplitude)` and turned the negative into `NaN`, collapsing the whole distribution to `Empty` under `--coalescent`; it now maps non-positive amplitudes to `+inf` cost (zero mass) via `Plain::is_defined`, so the negatives are tolerated. Any future consumer that assumes strictly non-negative `Function` amplitudes (a logarithm, a positivity assertion) would need the same guard.

## Options

- **Tolerate at consumers (current):** each operation that requires non-negative amplitudes treats `amplitude <= 0` as zero mass, matching `Plain::is_defined`. Robust but the invariant is enforced per consumer.
- **Clamp at the source:** clamp `multiply_function_function` output (and any other resampling producer) to `>= 0`, so `Function` amplitudes are non-negative by construction. Centralizes the invariant but changes the numeric output of every function-function product and is an unapproved numerical-behavior change.

## Related issues

- [M-timetree-backward-pass-plain-space-underflow.md](M-timetree-backward-pass-plain-space-underflow.md): the same product path underflows tail mass to zero; a separate scientific-accuracy concern.
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md): whether `NaN` amplitudes are rejected at construction; this issue removes one source of downstream `NaN`.
