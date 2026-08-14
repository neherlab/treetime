# Soft-tail boundary law evaluates in plain probability under negative-log storage

Timetree time distributions store negative-log ordinates $y = -\ln p$, where $p$ is the probability density [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81). The soft-tail boundary law `struct SoftTailLaw` still operates in plain probability, so it is on the wrong axis for those distributions.

## Mechanism

`fn SoftTailLaw::fit()` reads a grid ordinate and applies $-\ln$ to it, which treats the stored value as a plain probability [`packages/treetime-grid/src/soft_tail_law.rs#L50-L85`](../../packages/treetime-grid/src/soft_tail_law.rs#L50-L85). Its own doc comment states the precondition: the grid must hold plain-probability values [`packages/treetime-grid/src/soft_tail_law.rs#L45-L47`](../../packages/treetime-grid/src/soft_tail_law.rs#L45-L47). `fn SoftTailLaw::eval()` returns $p_\text{edge}\,e^{-k\,(t - t_\text{edge})}$ [`packages/treetime-grid/src/soft_tail_law.rs#L92-L94`](../../packages/treetime-grid/src/soft_tail_law.rs#L92-L94), where $p_\text{edge}$ is the edge probability, $t_\text{edge}$ the edge coordinate, and $k$ the neg-log slope.

The grid extrapolation seam passes the stored edge ordinate into `eval` as if it were $p_\text{edge}$, then stores the plain return value back as a neg-log ordinate [`packages/treetime-grid/src/grid_fn.rs#L378-L385`](../../packages/treetime-grid/src/grid_fn.rs#L378-L385). Under negative-log storage the edge ordinate is $-\ln p_\text{edge}$, not $p_\text{edge}$, so the evaluated tail is wrong.

The correct neg-log soft tail is a straight line in $t$:

$$y(t) = y_\text{edge} + k\,(t - t_\text{edge})$$

where $y_\text{edge} = -\ln p_\text{edge}$ is the stored edge ordinate. The sibling law `struct HardApproachLaw` was already migrated to this form: it fits and evaluates directly on the stored ordinate with no plain conversion [`packages/treetime-grid/src/hard_approach_law.rs#L57-L92`](../../packages/treetime-grid/src/hard_approach_law.rs#L57-L92). The soft law must adopt the same axis.

## Impact

The defect is latent. No production path fits a `SoftTailLaw` yet: the timetree passes attach only `Constant` and `Hard` boundaries, not `Linear`. The soft tail therefore blocks the log-space soft-tail step rather than producing a wrong result today. Any code that attaches a fitted `Linear` soft tail to a negative-log distribution would read a garbage extrapolation.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part A moves the axis to negative-log; Part B classifies `Linear` as the soft-boundary law.
- [N-timetree-inference-omits-integrable-soft-tail-boundaries.md](N-timetree-inference-omits-integrable-soft-tail-boundaries.md): the inference integration that this fix unblocks.
