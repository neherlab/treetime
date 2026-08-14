# Backward child fold multiplies per child instead of summing on a common grid

The timetree backward pass folds child time messages by multiplying and normalizing once per child. Under negative-log storage the proposal replaces this with a single common grid: choose one working grid for the node, resample each child message onto it once, sum, then re-window once at the end.

## Mechanism

The child fold calls `distribution_multiplication(&current, parent_message)?.normalize()` for every child [`packages/treetime/src/timetree/inference/backward_pass.rs#L82-L88`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L82-L88), and repeats the same multiplication and normalization for the date constraint [`packages/treetime/src/timetree/inference/backward_pass.rs#L104-L108`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L104-L108). Each `distribution_multiplication` builds its own grid over the support intersection and interpolates both operands onto it, and each `normalize()` runs between steps. A node with many children is therefore resampled and normalized once per child, in a fixed left-to-right order.

Under `NegLog` multiplication is addition, which is exact and associative with no intermediate normalization. Repeated per-child resampling adds interpolation drift and forces a resampling order.

## Required behavior

- Choose the node's working grid once: intersect the hard bounds, choose the soft extent generously, take spacing from the finest operand.
- Resample each child message and the date constraint onto that grid exactly once.
- Sum the ordinates in any order.
- Re-window and resample once, at the end.

Each message is then resampled exactly once regardless of fan-out, which removes the per-step drift and the order dependence and makes a balanced multiplication tree unnecessary.

## Impact

No demonstrated wrong result today; this is accuracy and structure, not a crash. It is proposal step 4 and the precondition for the adaptive mass-bounded domain, which re-windows once per node ([M-timetree-adaptive-mass-bounded-domain.md](M-timetree-adaptive-mass-bounded-domain.md)).

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part D, "Backward pass restructuring".
- [M-timetree-backward-pass-plain-space-underflow.md](M-timetree-backward-pass-plain-space-underflow.md): the plain-space underflow that motivated per-step normalization, closed by the negative-log switch.
