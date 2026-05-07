# Proposal: configurable convergence norm for GTR inference

## Summary

Add an option to select L2 or RMS convergence norm for site-specific GTR inference. Default L2 for v0 parity.

## Motivation

The iterative GTR inference uses L2 norm of the full pi difference array as convergence criterion. L2 norm scales with `sqrt(n_states * seq_len)`, so the per-element convergence threshold varies with problem size:

| seq_len | n_states | Per-element threshold (dp=1e-5) |
| ------: | -------: | ------------------------------: |
|       2 |        4 |                          3.5e-6 |
|     100 |        4 |                          5.0e-7 |
|   10000 |        4 |                          5.0e-8 |

For short sequences, the solver may stop with looser per-element accuracy. For long sequences, it over-converges and wastes iterations.

RMS norm (`sqrt(mean(delta^2))`) makes the threshold dimension-independent.

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

- For nucleotides (n_states=4, typical seq_len 100-30000): the L2 thresholds (5e-7 to 5e-8) are adequate. Practical impact is negligible.
- For amino acids (n_states=20) with long alignments: RMS would save iterations without loosening per-element accuracy.
- v0 uses L2. Changing the default would diverge from v0 convergence behavior.
