# Timetree passes attach a non-integrable Constant tail instead of an integrable Linear soft tail

The negative-log timetree passes attach a flat `Constant` tail on the soft side of every message. A `Constant` tail is a flat neg-log line, so its probability integral over the half-line diverges. The proposal retires `Constant` from the inference passes and attaches the integrable `Linear` soft tail instead.

## Mechanism

The passes set `Constant` on the unbounded side and `Hard(None)` on the bounded side after each message is built [`packages/treetime/src/timetree/inference/backward_pass.rs#L138-L140`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L138-L140) [`packages/treetime/src/timetree/inference/forward_pass.rs#L230-L231`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L230-L231) [`packages/treetime/src/timetree/inference/forward_pass.rs#L257-L258`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L257-L258). No pass produces a `Linear(Some(..))` tail, so `SoftTailLaw::fit` is never called in production and the whole `Linear` path is unused there.

`SoftTailLaw::fit` and `SoftTailLaw::eval` already operate on the negative-log stored-ordinate axis, and the `Linear` multiplication composition already exists and is correct [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L331-L340`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L331-L340). The remaining work is to fit and attach a `Linear` tail in the passes and to retire `Constant` from them.

## Required behavior

- Fit and attach a `Linear` soft tail on the soft side of each message, keeping the hard side hard.
- Retire `Constant` from the inference passes; keep it only for genuinely uninformative edges outside inference.
- The attached tail must survive normalization, which depends on the boundary law surviving regridding.

## Impact

No demonstrated runtime effect today: the confidence path that a `Constant` tail would corrupt is disabled under negative-log storage, and no `Linear` soft tail is produced yet. The gap is forward work toward the log-space soft-tail step, and it is a precondition for the mass-based domain, whose quantile integral needs an integrable tail.

## Related

- [M-grid-from-grid-array-discards-boundary-law.md](M-grid-from-grid-array-discards-boundary-law.md): prerequisite so the fitted tail survives normalization.
- [N-distribution-convolution-tail-not-refit-as-law.md](N-distribution-convolution-tail-not-refit-as-law.md): the convolution must produce a refit tail for the passes to consume.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part B retires `Constant`; design axis 4.
