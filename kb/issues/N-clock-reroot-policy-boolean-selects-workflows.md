# Clock reroot boolean selects incompatible workflows

`fn estimate_clock_model_with_reroot_policy()` [packages/treetime/src/clock/clock_regression.rs#L23-L28](../../packages/treetime/src/clock/clock_regression.rs#L23-L28) reads the field `keep_root: bool` of `struct ClockFitParams` [packages/treetime/src/clock/clock_regression.rs#L124-L130](../../packages/treetime/src/clock/clock_regression.rs#L124-L130) to select between backward regression alone and rerooting followed by backward regression.

`reroot_params` remains mandatory when `keep_root` makes it irrelevant. Call sites set the named field, but some still pass a literal, for example `keep_root: false` in [packages/treetime/src/timetree/optimization/reroot.rs#L68](../../packages/treetime/src/timetree/optimization/reroot.rs#L68), so the selected scientific workflow and the ignored parameters are visible only through the value of that field.

## Workflow difference

- Keep-root runs the backward regression on the input topology and returns no reroot result
- Optimize-root first calls `fn reroot_at_best_root()`, which searches for a root, mutates the topology, and returns reroot diagnostics; then it runs the same backward regression on the new topology [packages/treetime/src/clock/clock_regression.rs#L46-L78](../../packages/treetime/src/clock/clock_regression.rs#L46-L78)

The boolean therefore controls mutation, algorithm selection, and return-state availability. It also permits invalid request shapes: reroot configuration can be supplied while rerooting is disabled, and a caller can accidentally invert a `true` or `false` literal.

Represent root handling as an exhaustive policy. The keep-root state carries no reroot configuration; the optimize-root state carries the parameters required for topology mutation. Existing regression and reroot behavior must remain unchanged.

The result should likewise identify whether topology was retained or optimized instead of relying on an unrelated `Option` whose absence must be interpreted alongside the input flag.

## Validation

- Parameterized tests cover both policy variants and prove their allowed configuration differs by construction.
- Golden-master clock and timetree outputs preserve regression estimates, selected roots, and topology.
- A topology mutation assertion proves keep-root cannot invoke the root search.
- Call sites contain no boolean literals selecting root policy.

## Related issues

- [N-core-units-of-measurement-not-tracked.md](N-core-units-of-measurement-not-tracked.md)
