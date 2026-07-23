# GTR model state can be mutated into inconsistent configurations

`GTR::new()` normalizes rates and frequencies and derives `average_rate`, eigenvalues, eigenvectors, and the inverse eigenvector matrix. Every input and derived field remains publicly writable afterward [`packages/treetime/src/gtr/gtr.rs#L170-L201`](../../packages/treetime/src/gtr/gtr.rs#L170-L201).

Refinement writes `mu` directly through `gtr_mut()` [`packages/treetime/src/gtr/refinement.rs#L69`](../../packages/treetime/src/gtr/refinement.rs#L69) [`packages/treetime/src/gtr/refinement.rs#L148`](../../packages/treetime/src/gtr/refinement.rs#L148). `set_site_rates()` accepts any vector without checking length, finiteness, or positivity before propagation broadcasts it against sequence profiles [`packages/treetime/src/gtr/gtr.rs#L317`](../../packages/treetime/src/gtr/gtr.rs#L317).

## Invariants at risk

- `W` and `pi` are normalized together during construction.
- `average_rate`, `eigvals`, `v`, and `v_inv` are derived from the normalized rate matrix.
- `site_rates`, when present, must match sequence length and contain valid finite rates.
- Rate updates must preserve the coordinate convention used by propagation and branch optimization.

Public mutation allows a caller to change one member without updating the others. `expQt()` can then combine eigenvectors and eigenvalues derived from an older matrix. A site-rate length mismatch reaches ndarray broadcasting during `evolve()` instead of being rejected at construction.

`GTRSiteSpecific::new()` already validates rate length and rejects negative rates, so equivalent scientific validity currently depends on which model type is used.

## Required invariant

Callers must not be able to change normalized inputs or derived eigendecomposition state independently. Rate and site-rate updates need operations that validate their domain and either recompute dependent state or document why the derived state remains valid.

## Required design

- Make constructor inputs and derived state private.
- Expose read-only accessors for model inspection.
- Parse site rates into a validated value whose length is tied to the relevant sequence contract.
- Replace direct `mu` writes with an operation that states whether changing `mu` rescales branch coordinates or requires recomputation.
- Keep serialization through a validated constructor so deserialization cannot bypass invariants.

## Validation

- Reject non-finite, negative, and wrong-length site rates with exact contextual errors.
- Property tests verify normalized frequencies, symmetric exchangeabilities, eigendecomposition consistency, and stochastic transition matrices after every supported update.
- Golden-master comparisons cover refinement paths that currently mutate `mu`.
- Invalid intermediate model states are unrepresentable through the public API.

## Related issues

- [M-core-units-of-measurement-not-tracked.md](M-core-units-of-measurement-not-tracked.md)
- [M-gtr-per-site-rate-variation.md](M-gtr-per-site-rate-variation.md)
- [N-gtr-site-specific-partition-integration.md](N-gtr-site-specific-partition-integration.md)
