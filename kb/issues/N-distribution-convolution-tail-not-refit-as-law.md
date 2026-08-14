# Convolution extends tail samples but attaches no refit boundary law

The function-function convolution round-trips through a plain-space FFT and reconstructs decaying tail samples in negative-log space, but it attaches no fitted `SoftTailLaw` or `HardApproachLaw`. The proposal requires the convolved message to carry a refit law, because the FFT cannot produce the far tail and only a log-space fit reproduces it.

## Mechanism

The convolution normalizes each operand to plain space, runs the FFT, then reconstructs the tails by extrapolating the outermost trusted grid samples along a two-point secant in negative-log space [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L261-L319`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L261-L319). It extends grid ordinates only; it fits no law object and calls no `SoftTailLaw::fit` or `HardApproachLaw::fit`. The result therefore carries extended samples but an empty boundary law.

The caller then overwrites the tails with fixed classes and no law: the backward message gets `Constant` on the left and `Hard(None)` on the right [`packages/treetime/src/timetree/inference/backward_pass.rs#L138-L140`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L138-L140). Whatever the convolution reconstructed is discarded past the grid edge.

## Required behavior

- Fit the soft-tail law (and the hard-approach law where the edge is hard) in negative-log space from the outermost trusted points of the FFT output, below the roundoff floor.
- Attach the fitted law to the result so multiplication and re-windowing propagate it in closed form.
- Reproduce the reference guard: extend a tail only where it decays away from support; otherwise clamp.

## Impact

No demonstrated wrong result today: the confidence path that a non-integrable tail would corrupt is disabled under negative-log storage, and no fitted soft tail is produced yet. The gap is forward work. Its importance is that the decisive value in a disjoint-support product is far below peak, inside the region the FFT cannot produce but a log-linear fit reproduces exactly.

## Related

- [N-timetree-passes-omit-integrable-linear-soft-tail.md](N-timetree-passes-omit-integrable-linear-soft-tail.md): the passes must then consume the refit tail instead of overwriting it with `Constant`.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): "Accuracy constraint on convolution".
