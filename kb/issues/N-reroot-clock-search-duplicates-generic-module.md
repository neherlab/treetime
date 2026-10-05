# Clock rerooting duplicates the generic root-search implementation

The reusable generic root-search module implements pluggable `RootStats`, edge cost evaluation, split optimization, and topology orchestration under [`packages/treetime/src/reroot`](../../packages/treetime/src/reroot). Clock rerooting still maintains a parallel implementation under [`packages/treetime/src/clock/find_best_root`](../../packages/treetime/src/clock/find_best_root), so fixes and invariants can drift between two implementations of the same search mechanics.

Both modules also export `FindRootResult`, `find_best_root()`, and `find_best_split()` names for different result contracts. Imports and diagnostics therefore conceal whether a value carries clock regression statistics or generic objective statistics.

The migration must preserve payload statistics across the module boundary and use a mechanical-equivalence oracle for every clock and timetree caller. Optimizer defaults remain outside that oracle because changing them alters numerical behavior.

## Potential solutions

- O1. Implement `RootStats` directly for `ClockSet` and use the generic engine without an adapter.
- O2. Introduce a dedicated `ClockRootStats` value that owns only the sufficient statistics required by the generic engine. This makes search inputs explicit but duplicates part of `ClockSet`'s representation.

## Recommendation

Use O1. Implement `RootStats` directly for `ClockSet`, route clock rerooting through the generic search, preserve the existing clock objective and explicitly supplied optimizer parameters, and delete the duplicate clock-only search module after all callers migrate. Objective/default changes and tip-name policy remain separate issues.

## Fix (O1)

Route the existing clock reroot objective through the generic `reroot` infrastructure without changing numerical behavior, defaults, or name-resolution policy:

- Implement `trait RootStats` [packages/treetime/src/reroot/traits.rs#L3](../../packages/treetime/src/reroot/traits.rs#L3) directly for `struct ClockSet` [packages/treetime/src/clock/clock_set.rs#L9](../../packages/treetime/src/clock/clock_set.rs#L9) through `leaf_contribution_to_parent()`, `propagate_averages()`, and `chisq()` (`clock_set.rs#L41`, `#L55`, `#L122`)
- Replace the clock-only search and edge-cost calls with the corresponding generic `reroot` search and orchestration APIs
- Extract the existing per-edge `ClockSet` statistics into the map that the generic search reads
- Preserve every caller-supplied and default `BranchPointOptimizationParams` value exactly
- Preserve the fixed-zero-rate MinDev objective (`RootObjective::FixedRate(0.0)` [packages/treetime/src/clock/reroot.rs#L194-L205](../../packages/treetime/src/clock/reroot.rs#L194-L205)) and the current tip-name resolution behavior
- Keep the clock-specific fixup that swaps the `clock_to_parent` and `clock_to_child` statistics on inverted edges, in `fn apply_reroot()` [packages/treetime/src/clock/reroot.rs#L347-L371](../../packages/treetime/src/clock/reroot.rs#L347-L371)
- Update imports and delete `packages/treetime/src/clock/find_best_root/` after no callers remain

## Validation

- For explicit grid, Brent, and golden-section parameters, compare the selected edge, split fraction, score, clock model, and output before and after the migration
- Cover least-squares, min-dev, and tip reroot modes on deterministic fixtures
- All existing clock and timetree reroot tests pass without changes to expected values
- Unit-test the `RootStats` methods against the corresponding `ClockSet` methods over the same inputs

## Related issues

- [N-reroot-duplicated-tip-name-resolution.md](N-reroot-duplicated-tip-name-resolution.md)
- [N-reroot-split-optimizer-default-diverges-from-v0.md](N-reroot-split-optimizer-default-diverges-from-v0.md)

## Related errata

- [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)
