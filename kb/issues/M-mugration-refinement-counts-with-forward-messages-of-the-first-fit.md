# Mugration refinement counts with the forward messages of the first fitted model

`fn refine_gtr()` [packages/treetime/src/gtr/refinement.rs#L17-L83](../../packages/treetime/src/gtr/refinement.rs#L17-L83) runs one full marginal update under the first fitted model and passes its forward messages and root state to `count_transitions` in every later iteration ([packages/treetime/src/gtr/refinement.rs#L37-L68](../../packages/treetime/src/gtr/refinement.rs#L37-L68)). The loop replaces the backward messages after each rate step but does not recompute the forward messages, so from the second fit on the counts combine backward messages and transition matrices of the current model with forward messages of the first fitted model and rate.

This is the v0 order: v0's `infer_gtr(marginal=True)` reuses its existing marginal reconstruction ([packages/legacy/treetime/treetime/treeanc.py#L1544](../../packages/legacy/treetime/treetime/treeanc.py#L1544)), and the mugration golden masters match v0 within 1e-7 with it ([kb/algo/mugration.md](../algo/mugration.md#golden-master-tolerance)). Counts from messages of different models are not an expectation step of EM, so the log likelihood is not guaranteed to increase from one iteration to the next.

> [!IMPORTANT]
> **Decision required.** Options:
>
> - Keep the v0 order and close this issue
> - Run a full marginal update under the current model before each fit, as [kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md) proposes, and record the divergence from v0. The golden masters then stop matching v0 within 1e-7

## Validation

- With the full update per iteration: the log likelihood does not decrease from one iteration to the next on the golden-master datasets
