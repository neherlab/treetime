# Profile resolution boolean combines deterministic and stochastic operations

`fn resolve_profile()` accepts `sample: bool` and an RNG for both deterministic argmax resolution and posterior sampling [`packages/treetime/src/ancestral/sample.rs#L49`](../../packages/treetime/src/ancestral/sample.rs#L49). The deterministic path does not use the RNG, while the stochastic path consumes mutable RNG state.

## Current contract

- `sample == true` delegates to `sample_from_profile()` and advances the caller's RNG.
- `sample == false` delegates to deterministic `argmax_first()` but still requires a mutable RNG argument.
- Call sites pass a bare boolean, so the operation and its reproducibility implications are invisible without reading the callee.

The two branches have different observable contracts: one is deterministic and pure with respect to RNG state; the other is stochastic and stateful. Combining them forces every caller to satisfy the union of both interfaces.

## Required boundary

Expose separate deterministic-resolution and posterior-sampling operations. Select the operation where reconstruction policy is known, and require an RNG only for sampling. Preserve deterministic first-maximum tie breaking and the existing sampling distribution.

## Validation

- Unit cases prove deterministic tie handling without constructing an RNG.
- Seeded sampling tests verify the same categorical distribution and reproducibility contract.
- Call-site tests make sampling policy explicit for dense, sparse, and tip-state reconstruction paths.
