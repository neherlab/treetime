# Empty-from-disjoint-hard invariant is enforced in multiplication only

The proposal makes `Distribution::Empty` reachable only from genuinely disjoint hard domains, as a checked error rather than a silent collapse. Multiplication enforces this, but convolution, division, and the point/range/formula multiplication arms still return `Empty` silently.

## Mechanism

Multiplication routes an empty result through `multiplication_empty_result`, which returns `Empty` when `hard_domains_disjoint` holds and raises `make_internal_error!` otherwise [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L406-L456`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L406-L456). This guard is reached only from the `Disjoint` arm of the three function-producing multiplies. The remaining multiply arms (point/point, point/range, range/range, formula variants) and the `Empty` propagation arm still return `Distribution::empty()` unguarded.

Convolution returns `Empty` at four sites with no guard [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L169`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L169) [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L182`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L182) [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L198`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L198) [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L206`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L206). Division returns `Empty` on disjoint support with no guard [`packages/treetime-distribution/src/distribution_ops/divide.rs#L93`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L93) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L104`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L104) [`packages/treetime-distribution/src/distribution_ops/divide.rs#L126`](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L126).

## Required behavior

- Route every `Empty` producer through the same check: `Empty` is permitted only when the operands' hard domains are genuinely disjoint; any other empty result is an internal error.
- Reuse `hard_domains_disjoint` so the rule is single-sourced across multiplication, convolution, and division.

## Impact

Silent `Empty` from numerical collapse was the original motivating defect: an empty message poisons every ancestor to the root and prints `log_lh_pos = NaN` with no warning. The guard closes that path for the function multiplies but leaves the other operators able to reintroduce it. This is proposal step 8.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part B "Invariant" and design axis 7.
