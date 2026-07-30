# Timetree silently drops node dates when time distribution is empty

When an internal node's time distribution is `Empty` (grid length 0), `fn Distribution.likely_time()` returns `None`, and `fn set_likely_time()` at [packages/treetime/src/timetree/inference/forward_pass.rs#L58-L69](../../packages/treetime/src/timetree/inference/forward_pass.rs#L58-L69) silently skips setting the node's time. The output tree contains no date fields for the affected nodes -- no error, no warning.

## Impact

An empty time distribution indicates that the inference constraints at a node are irreconcilable (e.g. disjoint child messages). The silent skip produces an output tree that looks structurally valid but is missing all date information for internal nodes. Users have no signal that the inference failed.

## Suggested fix

Emit a warning or error when `likely_time()` returns `None` for an internal node after the forward pass. At minimum, a log message naming the affected node. A stricter version would return an error to the caller so the pipeline can report the failure.

## Related

The most common cause of empty distributions -- disjoint backward-message multiplication under `--keep-root` -- was fixed by honoring `Constant` tails in multiplication. This issue is about the safety net: even with the tail fix, other scenarios (severely conflicting date constraints, numerical edge cases) could produce empty distributions, and the silent skip would hide them.
