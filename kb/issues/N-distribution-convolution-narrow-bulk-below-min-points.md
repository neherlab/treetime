# Convolution rejects a valid result when the trusted bulk resamples below two points

The function-function convolution crops its FFT output to the trusted bulk, then resamples to the coarser operand spacing. A very narrow bulk can resample to fewer than two points, and the routine returns an error instead of the narrow result.

## Mechanism

After the FFT and tail reconstruction, the routine resamples the result to the coarser of the two operand spacings, then errors when the point count is below two [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L223-L229`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L223-L229). The trusted bulk is the region above the roundoff floor [`packages/treetime-distribution/src/distribution_ops/convolve.rs#L270-L302`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L270-L302). When two narrow operands produce a bulk only a few fine cells wide, the coarse resample can collapse it below two points. A single-cell result before the resample is already handled as a point mass, so this path is specific to the resample step.

## Impact

The defect is latent. No test drives a convolution narrow enough to reach it. The crop to the trusted bulk makes the path marginally more reachable than before, because the kept range can be shorter than the raw convolution. A narrow but valid convolution would fail rather than return a coarse two-point or point-mass result.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part D makes the grid a mass-bounded domain with a resolution floor, which is the context for choosing the minimum point count after a resample.
