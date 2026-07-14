# Clock rerooting duplicates the generic root-search implementation

The reusable generic root-search module implements pluggable `RootStats`, edge cost evaluation, split optimization, and topology orchestration under [`packages/treetime/src/reroot`](../../packages/treetime/src/reroot). Clock rerooting still maintains a parallel implementation under [`packages/treetime/src/clock/find_best_root`](../../packages/treetime/src/clock/find_best_root), so fixes and invariants can drift between two implementations of the same search mechanics.

## Potential solutions

- O1. Implement `RootStats` directly for `ClockSet` and use the generic engine without an adapter.
- O2. Introduce a dedicated `ClockRootStats` value that owns only the sufficient statistics required by the generic engine. This makes search inputs explicit but duplicates part of `ClockSet`'s representation.

## Recommendation

Use O1. Implement `RootStats` directly for `ClockSet`, route clock rerooting through the generic search, preserve the existing clock objective and explicitly supplied optimizer parameters, and delete the duplicate clock-only search module after all callers migrate. Objective/default changes and tip-name policy remain separate issues.

## Related issues

- [M-clock-mindev-wrong-objective.md](M-clock-mindev-wrong-objective.md)
- [N-reroot-duplicated-tip-name-resolution.md](N-reroot-duplicated-tip-name-resolution.md)
- [N-reroot-split-optimizer-default-diverges-from-v0.md](N-reroot-split-optimizer-default-diverges-from-v0.md)
