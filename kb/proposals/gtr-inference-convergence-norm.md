# Proposal: configurable convergence norm for GTR inference

## Summary

Add an option to select L2 or RMS convergence norm for site-specific GTR inference. Default L2 for v0 parity.

## Motivation

The iterative GTR inference uses the L2 norm of the full $\pi$-difference array as its convergence criterion. For $N=\text{n\_states}\times\text{seq\_len}$, the equivalent RMS cutoff is $\mathit{dp}/\sqrt{N}$:

| seq_len | n_states | Equivalent RMS cutoff (`dp=1e-5`) |
| ------: | -------: | ------------------------------: |
|       2 |        4 |                          3.5e-6 |
|     100 |        4 |                          5.0e-7 |
|   10000 |        4 |                          5.0e-8 |

This conversion describes aggregate RMS change; it is not a bound on every element. A fixed RMS threshold with the same numeric `dp` would be looser than the current L2 criterion by a factor of $\sqrt{N}$.

RMS norm (`sqrt(mean(delta^2))`) makes the aggregate cutoff dimension-independent, but selecting its default value is a separate numerical decision.

## Design

```rust
#[derive(Clone, Debug, SmartDefault)]
pub enum ConvergenceNorm {
  #[default]
  L2,
  Rms,
}
```

Add to `InferGtrSiteSpecificOptions`. Default `L2` preserves v0 behavior. The change in the solver is: divide `dist` by `sqrt(n_elements)` when `Rms`.

## Impact

- No current benchmark establishes the effect on iterations, runtime, inferred parameters, or likelihood.
- Using RMS with the same numeric cutoff loosens the aggregate stopping criterion relative to L2; parity therefore requires L2 as the default.
- v0 uses L2. Changing the default would diverge from v0 convergence behavior.
