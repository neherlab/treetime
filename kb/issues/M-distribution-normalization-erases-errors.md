# Distribution normalization erases formula and grid errors

Negative-log normalization converts formula-evaluation and grid-construction failures into `Distribution::Empty`. Callers therefore cannot distinguish an explicitly empty domain value from a failed likelihood calculation.

## Evidence

`Distribution<NegLog>::to_plain_normalized()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L378](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L378) uses `map_or(Distribution::Empty, ...)` when formula discretization fails. `fn neglog_function_to_plain_normalized()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L437](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L437) also maps `DistributionFunction::from_start_dx_values()` errors to `Distribution::Empty`.

These conversions feed the timetree backward pass in [packages/treetime/src/timetree/inference/backward_pass.rs#L58](../../packages/treetime/src/timetree/inference/backward_pass.rs#L58), where the empty value can be stored as a successful node-time distribution.

## Required contract

- Normalization returns `Result<Distribution<Plain>, Report>`.
- `Distribution::Empty` maps to `Ok(Distribution::Empty)` because emptiness is a domain value.
- Formula evaluation and grid construction propagate their original error with distribution context.
- Callers propagate the error and never store a failed normalization as an empty time distribution.

## Potential solutions

- O1. Return `Result<Distribution<Plain>, Report>` and reserve `Distribution::Empty` for domain emptiness.
- O2. Add an error-bearing distribution variant. This keeps a non-fallible signature but allows failed computations to travel as values through likelihood code.

## Recommendation

Use O1: make normalization fallible and reserve `Distribution::Empty` for explicit domain emptiness. The independent policy for non-finite sampled values remains tracked separately.

## Required properties

For finite negative-log samples $\ell_i$, max-normalized weights are

$$w_i = e^{-(\ell_i-\ell_{\min})}$$

where $\ell_{\min}=\min_j\ell_j$. The output is finite and non-negative, satisfies $\max_i w_i=1$, preserves $w_i/w_k=e^{-(\ell_i-\ell_k)}$ for representable nonzero weights, and is invariant under adding a common finite offset to every $\ell_i$.

## Related issues

- [N-distribution-formula-silent-discretization.md](N-distribution-formula-silent-discretization.md)
- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md)
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md)
