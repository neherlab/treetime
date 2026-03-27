# Documentation

## Directories

- [`docs/algorithms/`](algorithms/) - Design documents specifying intended v1 capabilities. Features described here are work items even if v0 does not implement them.
- [`docs/port-algo-inventory/`](port-algo-inventory/_index.md) - Scientific background, algorithm catalog, v0 locations. [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) has v0 algorithm detail for porting.
- [`docs/port-feature-inventory/`](port-feature-inventory/_index.md) - v0/v1 feature parity checklist. `[x]` implemented, `[/]` partial, `[ ]` missing.
- [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) - Deliberate v1 deviations from v0, or intentional removals, with rationale.
- [`docs/port-known-issues/`](port-known-issues/_index.md) - All unfinished work items by severity (C/H/M/N prefix): bugs, missing v0 features, missing design-doc features, stubs, dead flags. Includes [`docs/port-known-issues/_test-matrix.md`](port-known-issues/_test-matrix.md) for systematic v0/v1 comparison.
- [`docs/port-proposals/`](port-proposals/_index.md) - New v1 features that neither v0 nor the design docs specify, pre-implementation.
- [`docs/port-test-inventory/`](port-test-inventory/_index.md) - Test coverage by domain.
- [`docs/port-v0-errata/`](port-v0-errata/_index.md) - Defects and oversights in v0 that v1 does not reproduce. ALWAYS check before matching v0 behavior.
- [`docs/reports/`](reports/) - Research reports on iterative tree refinement and optimization methods.
- [`docs/dev/`](dev/) - Development infrastructure documentation.
- [`docs/cli_overview.md`](cli_overview.md) - CLI command overview.

## Porting Ledger Taxonomy

Source code is ground truth. Ledger entries are guides - NEVER substitute for code verification.

Work items come from three sources: v0 code (porting gaps), `docs/algorithms/` design documents (specified features), and testing (discovered bugs). All three feed into `port-known-issues/`.

Every work item falls into exactly one category:

| Category                        | Directory                                                                                                                                     | Scope                                                                                                                                                                                                                                               |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Implemented, same behavior      | [`docs/port-algo-inventory/`](port-algo-inventory/_index.md), [`docs/port-feature-inventory/`](port-feature-inventory/_index.md)              | v0 or design-doc feature ported with equivalent behavior                                                                                                                                                                                            |
| Implemented, different behavior | [`docs/port-intentional-changes/`](port-intentional-changes/_index.md)                                                                        | v0 feature ported with deliberate divergence, or intentionally removed (one file per deviation)                                                                                                                                                     |
| Not yet done                    | [`docs/port-known-issues/`](port-known-issues/_index.md), [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) | Bugs, missing v0 features, missing design-doc features, stubs, dead flags. Known issues cover all types; [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md) adds algorithm-specific detail for unported algorithms |
| New in v1                       | [`docs/port-proposals/`](port-proposals/_index.md)                                                                                            | Feature that neither v0 nor the design docs specify. Moves to [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) after implementation                                                                                           |
| v0 defective, v1 correct        | [`docs/port-v0-errata/`](port-v0-errata/_index.md)                                                                                            | Defect or oversight in v0 that v1 does not reproduce (one file per erratum with evidence)                                                                                                                                                           |

### Decision rules

- v0 or design-doc feature working in v1 with same results -> [`docs/port-feature-inventory/`](port-feature-inventory/_index.md) `[x]`, algorithm in domain file
- v0 feature working in v1 with different results -> [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (one file with rationale)
- v0 feature intentionally removed from v1 -> [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (one file with rationale)
- v0 feature missing, stubbed, or broken in v1 -> [`docs/port-known-issues/`](port-known-issues/_index.md) (severity-prefixed file). If algorithmic, also in [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md)
- Design-doc feature not yet implemented -> [`docs/port-known-issues/`](port-known-issues/_index.md) (severity-prefixed file). If algorithmic, also in [`docs/port-algo-inventory/unimplemented.md`](port-algo-inventory/unimplemented.md)
- v1 feature that neither v0 nor design docs specify -> [`docs/port-proposals/`](port-proposals/_index.md) (pre-implementation), then [`docs/port-intentional-changes/`](port-intentional-changes/_index.md) (post-implementation)
- v0 behavior appears wrong, v1 does the correct thing -> [`docs/port-v0-errata/`](port-v0-errata/_index.md) (one file with evidence). Distinct from intentional changes (v1 design choices) and known issues (v1 defects)
- NEVER put v0 port targets or design-doc features in [`docs/port-proposals/`](port-proposals/_index.md)
