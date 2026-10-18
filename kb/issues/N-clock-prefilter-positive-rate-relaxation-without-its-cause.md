# Clock pre-filter relaxes the positive rate for a failure the exact root search removed

## Problem

The clock pre-filter searches the root without requiring a positive rate
([clock-prefilter-relaxed-positive-rate.md](../decisions/clock-prefilter-relaxed-positive-rate.md)), where v0 requires
one. The decision's rationale is that v1 rejected every root of `data/dengue/100` for a negative rate before outlier
filtering, so requiring a positive rate stopped the run.

That failure came from the old root search, which optimized the split only on the edges next to the best node. With
the exact search on every edge ([clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)),
`treetime clock` on `data/dengue/100` with `--clock-filter 0` (all tips in the fit, positive rate required) finds the
root at rate $6.88\cdot10^{-4}$; binaries from before the exact search stop with "Clock rate is negative for all root
positions" on the same command.

> [!IMPORTANT]
> **Decision required.** Restore v0's positive-rate requirement in the pre-filter, which the original reason no
> longer argues against, or keep the relaxation with a new rationale (for example, that the outlier screen should
> not stop a run whose final fit would succeed). Restoring it can change the pre-filter root and the outlier sets.

## Validation

With the requirement restored: the smoke matrix and the 240 validation sets with `--clock-filter` run without new
failures, and the outlier sets are compared against the current ones.
