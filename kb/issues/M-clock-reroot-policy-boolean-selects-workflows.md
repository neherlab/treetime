# Clock reroot boolean selects incompatible workflows

`fn estimate_clock_model_with_reroot_policy()` accepts `keep_root: bool` to select between backward regression alone and backward regression followed by forward regression and topology mutation [`packages/treetime/src/clock/clock_regression.rs#L238`](../../packages/treetime/src/clock/clock_regression.rs#L238).

Reroot and optimization parameters remain mandatory when `keep_root` makes them irrelevant. Call sites pass raw booleans, so the selected scientific workflow and ignored arguments are not visible at the call site.

## Workflow difference

- Keep-root performs backward regression and returns no reroot result.
- Optimize-root performs forward regression, searches for a root, mutates topology, and returns reroot diagnostics.

The boolean therefore controls mutation, algorithm selection, and return-state availability. It also permits invalid request shapes: reroot configuration can be supplied while rerooting is disabled, and a caller can accidentally invert a bare `true` or `false`.

Represent root handling as an exhaustive policy. The keep-root state carries no reroot configuration; the optimize-root state carries the parameters required for topology mutation. Existing regression and reroot behavior must remain unchanged.

The result should likewise identify whether topology was retained or optimized instead of relying on an unrelated `Option` whose absence must be interpreted alongside the input flag.

## Validation

- Parameterized tests cover both policy variants and prove their allowed configuration differs by construction.
- Golden-master clock and timetree outputs preserve regression estimates, selected roots, and topology.
- A topology mutation assertion proves keep-root cannot invoke forward search.
- Call sites contain no boolean literals selecting root policy.

## Related issues

- [N-reroot-clock-search-duplicates-generic-module.md](N-reroot-clock-search-duplicates-generic-module.md)
- [M-core-units-of-measurement-not-tracked.md](M-core-units-of-measurement-not-tracked.md)
