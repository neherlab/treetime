# Distribution support-boundary semantics are unresolved

`GridFn` is a generic interpolated function with constant extrapolation, while `Distribution::Function` is also used as a finite-support probability distribution. These contracts conflict: a bounded probability distribution is effectively zero outside its support in v0, but v1 interpolation repeats the nearest boundary value indefinitely.

## Current and reference behavior

[`GridFn::interp()`](../../packages/treetime-grid/src/grid_fn.rs) returns the first or last stored value outside the grid. Every `DistributionFunction::interp()` and resampling operation inherits that behavior.

V0 [`Distribution`](../../packages/legacy/treetime/treetime/distribution.py) stores neg-log probability and configures `interp1d(..., fill_value=BIG_NUMBER)`, which is effectively zero probability outside `[xmin, xmax]`. V0 does not encode convolution-tail behavior as generic distribution extrapolation: [`NodeInterpolator.convolve_fft()`](../../packages/legacy/treetime/treetime/node_interpolator.py) separately constructs and conditionally extrapolates convolution tails from their slopes.

Commit `542ac860c7cfa4bab6764aee1d1b3810a09eb54f` proposes pass-specific constant tails on `GridFn`. That proposal conflicts with finite-support multiplication/division, where operations are defined on the support intersection. It is also representation-dependent: numeric `0.0` means zero probability under `Plain`, but it is the multiplicative identity (probability one) under `NegLog`. A generic `Zero` tail is therefore not policy-correct.

## Independent decision axes

### A1. Abstraction contract

- **Generic function:** `GridFn` keeps mathematical constant extrapolation; probability support is enforced by `DistributionFunction` or distribution operations.
- **Support-aware grid:** `GridFn` stores explicit left/right boundary behavior, affecting every non-probability caller and every transformation such as resampling and argument negation.

### A2. Out-of-support value by representation

- **Probability zero:** return `0.0` for `Plain` and positive infinity for `NegLog`.
- **Boundary continuation:** repeat the boundary value, meaning the represented distribution has an unbounded flat tail.
- **Evaluation error:** reject out-of-support evaluation and require each operation to establish a valid domain first.

### A3. Tail ownership

- **Distribution-level:** a distribution declares its own tail semantics everywhere.
- **Operation-level:** multiplication/division use finite support, while convolution explicitly constructs justified tails as v0 does.
- **Inference-pass-level:** backward and forward timetree passes override tails. This makes the same distribution mean different things depending on the caller and requires a scientific justification for each side.

### A4. Default parity target

- **V0 parity:** finite-support distributions are zero outside support; convolution-tail extension remains an explicit convolution operation.
- **Approved v1 divergence:** define and validate a different policy, including its effect on posterior normalization, node times, and topology constraints.

These axes can be combined independently. No combination is approved, so there is no implementation ticket.

The reference-aligned candidate is A1 generic function, A2 representation-aware probability zero at the distribution boundary, A3 operation-owned convolution tails, and A4 v0 parity. This follows the project's default parity rule and keeps convolution extrapolation explicit. It remains an unapproved scientific/numerical choice.

## Related issues

- [M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md](M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
