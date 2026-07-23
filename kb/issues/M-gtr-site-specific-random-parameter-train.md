# Site-specific GTR random construction exposes a parameter train

`fn GTRSiteSpecific::random()` accepts six model-generation values plus an RNG as positional arguments [`packages/treetime/src/gtr/gtr_site_specific.rs#L198`](../../packages/treetime/src/gtr/gtr_site_specific.rs#L198). The same type already uses `GTRSiteSpecificParams` for ordinary construction [`packages/treetime/src/gtr/gtr_site_specific.rs#L12`](../../packages/treetime/src/gtr/gtr_site_specific.rs#L12).

The generation values form one coherent configuration and are difficult to verify by position at call sites. Represent random generation with a named parameter object and keep the RNG as the separate stateful collaborator.

## Parameter roles

The positional values define state count, sequence length, average rate, Dirichlet concentration for equilibrium profiles, Dirichlet concentration for exchangeabilities, and Gamma shape for site rates. Several share the same primitive type, so transposition compiles while changing the generated model distribution.

Ordinary construction already demonstrates a parameter-object boundary. Random-generation parameters differ in purpose and should use their own validated type rather than overloading ordinary model state.

## Required contract

Parse generation inputs into a complete immutable configuration that validates dimensions and positive finite distribution parameters. The random constructor consumes that configuration and a seeded RNG, returning a fully validated `GTRSiteSpecific`.

## Validation

- Each invalid parameter has an exact rejection case at the configuration boundary.
- Seeded generation remains reproducible for the same configuration and PRNG.
- Property tests verify generated frequencies normalize, rates satisfy the configured domain, and transition matrices retain GTR invariants.
- Call sites contain no positional run of model-generation scalars.

## Related issues

- [N-gtr-site-specific-partition-integration.md](N-gtr-site-specific-partition-integration.md)
