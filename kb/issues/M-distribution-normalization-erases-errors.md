# Distribution normalization erases formula and grid errors

Negative-log normalization converts formula-evaluation and grid-construction failures into `Distribution::Empty`. Callers therefore cannot distinguish an explicitly empty domain value from a failed likelihood calculation.

## Evidence

`Distribution<NegLog>::normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L177-L199](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L177-L199) uses `map_or(Distribution::Empty, ...)` when formula discretization fails. `fn neglog_function_normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207) also returns `Distribution::Empty` when the minimum reduction fails or is not finite.

The timetree backward pass stores the normalized value as the node-time distribution in [packages/treetime/src/timetree/inference/backward_pass.rs#L41-L44](../../packages/treetime/src/timetree/inference/backward_pass.rs#L41-L44), so an empty value from a failed normalization is stored as a successful result.

## Required contract

- Normalization returns `Result<Distribution<NegLog>, Report>`.
- `Distribution::Empty` maps to `Ok(Distribution::Empty)` because emptiness is a domain value.
- Formula evaluation and grid construction propagate their original error with distribution context.
- Callers propagate the error and never store a failed normalization as an empty time distribution.

## Potential solutions

- O1. Return `Result<Distribution<NegLog>, Report>` and reserve `Distribution::Empty` for domain emptiness.
- O2. Add an error-bearing distribution variant. This keeps a non-fallible signature but allows failed computations to travel as values through likelihood code.

## Recommendation

Use O1: make normalization fallible and reserve `Distribution::Empty` for explicit domain emptiness. The independent policy for non-finite sampled values remains tracked separately.

## Required properties

For finite negative-log samples $\ell_i$, the normalized samples are

$$\ell'_i = \ell_i-\ell_{\min}$$

where $\ell_{\min}=\min_j\ell_j$. The output is finite and non-negative, satisfies $\min_i \ell'_i=0$, preserves the differences $\ell'_i-\ell'_k=\ell_i-\ell_k$, and is invariant under adding a common finite offset to every $\ell_i$.

## Related issues

- [N-distribution-formula-silent-discretization.md](N-distribution-formula-silent-discretization.md)
- [N-error-suppression-unwrap-or-defaults.md](N-error-suppression-unwrap-or-defaults.md)
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md)
