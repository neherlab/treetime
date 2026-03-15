# Documentation

## Directories

- [`docs/port-algo-inventory/`](port-algo-inventory/_index.md) - Scientific background, algorithm catalog, v0 locations. [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) has v0 algorithm detail for porting.
- [`docs/port-feature-inventory/`](port-feature-inventory/_index.md) - v0/v1 feature parity checklist. `[x]` implemented, `[/]` partial, `[ ]` missing.
- [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) - Deliberate v1 deviations from v0, or intentional removals, with rationale.
- [`docs/port-known-issues/`](port-known-issues/_index.md) - Bugs and missing features by severity (C/H/M/N prefix). Includes [`docs/port-known-issues/_test-matrix.md`](port-known-issues/_test-matrix.md) for systematic v0/v1 comparison.
- [`docs/port-proposals/`](port-proposals/_index.md) - New v1 features not in v0, pre-implementation.
- [`docs/port-test-inventory/`](port-test-inventory/_index.md) - Test coverage by domain.
- [`docs/algorithms/`](algorithms/) - Algorithm design notes (volatile, working documents).
- [`docs/dev/`](dev/) - Development infrastructure documentation.
- [`docs/cli_overview.md`](cli_overview.md) - CLI command overview.

## Porting Ledger Taxonomy

Source code is ground truth. Ledger entries are guides - NEVER substitute for code verification.

Every v0/v1 difference falls into exactly one category:

| Category                        | Directory                                                                                                                                     | Scope                                                                                                                                                                                                    |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Implemented, same behavior      | [`docs/port-algo-inventory/`](port-algo-inventory/_index.md), [`docs/port-feature-inventory/`](port-feature-inventory/_index.md)              | v0 feature ported with equivalent behavior                                                                                                                                                               |
| Implemented, different behavior | [`docs/port-intentional-changes/`](port-intentional-changes/_index.md)                                                                        | v0 feature ported with deliberate divergence, or intentionally removed (one file per deviation)                                                                                                          |
| Not yet ported                  | [`docs/port-known-issues/`](port-known-issues/_index.md), [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) | v0 feature missing from v1: bugs, stubs, dead flags. Known issues cover all types; [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) adds algorithm-specific v0 detail |
| New in v1                       | [`docs/port-proposals/`](port-proposals/_index.md)                                                                                            | Feature v0 does not have. Moves to [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) after implementation                                                                           |

### Decision rules

- v0 feature working in v1 with same results -> [`docs/port-feature-inventory/`](port-feature-inventory/_index.md) `[x]`, algorithm in domain file
- v0 feature working in v1 with different results -> [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (one file with rationale)
- v0 feature intentionally removed from v1 -> [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (one file with rationale)
- v0 feature missing, stubbed, or broken in v1 -> [`docs/port-known-issues/`](port-known-issues/_index.md) (severity-prefixed file). If algorithmic, also in [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md)
- v1 feature that v0 does not have -> [`docs/port-proposals/`](port-proposals/_index.md) (pre-implementation), then [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (post-implementation)
- NEVER put v0 port targets in [`docs/port-proposals/`](port-proposals/_index.md)
