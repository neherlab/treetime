# Timetree inference omits integrable soft-tail boundaries and discards convolution tails

The negative-log timetree passes attach flat and hard boundaries only. They never attach the integrable `Linear` soft-tail law, and they overwrite the tail that the convolution reconstructs. The log-space proposal requires the opposite on both points.

## Mechanism

Two mechanisms combine.

- The passes set `Constant` and `Hard(None)` boundaries after each message is built [`packages/treetime/src/timetree/inference/backward_pass.rs#L139-L140`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L139-L140) [`packages/treetime/src/timetree/inference/forward_pass.rs#L230-L231`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L230-L231) [`packages/treetime/src/timetree/inference/forward_pass.rs#L257-L258`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L257-L258). A `Constant` tail is a flat neg-log line, so its probability integral over the half-line diverges.
- The FFT convolution reconstructs a log-linear decaying tail from the outermost trusted points and bakes it as negative-log samples [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L261-L319`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L261-L319). It attaches no fitted `SoftTailLaw`, and the boundary the passes then set replaces whatever the convolution produced.

## Impact

The proposal retires `Constant` from the inference passes, because a non-integrable tail corrupts the mass integral behind the quantile and confidence-interval path. It replaces `Constant` with `Linear`, which is integrable in closed form. It also requires the convolution result to carry a refit soft-tail law that the passes consume, so that no tail mass leaks across operations.

The gap has no demonstrated runtime effect today. The confidence path that a `Constant` tail would corrupt is disabled under negative-log storage, and no `Linear` soft tail is produced yet. The gap is forward work toward the log-space soft-tail step.

The `Linear` multiplication composition already exists and is correct, so the boundary law is ready to carry once it is produced [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L331-L340`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L331-L340).

## Related

- [M-distribution-soft-tail-law-plain-space-under-neglog.md](M-distribution-soft-tail-law-plain-space-under-neglog.md): the soft-tail law must first evaluate on the negative-log axis. This issue is blocked on it.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part B retires `Constant`; the convolution accuracy constraint requires a refit tail from the trusted points.
