# Marginal timetree point estimates use an unapproved topology projection

Marginal inference computes each node's posterior independently, so the raw posterior modes are not jointly constrained by `child_time >= parent_time` in calendar coordinates. Current v1 then projects the committed point estimate of each non-leaf internal node onto that constraint. The projection prevents negative committed durations between internal nodes but does not make the marginal posteriors jointly consistent.

## Current and reference behavior

After refining a node's marginal distribution, [`set_likely_time()`](../../packages/treetime/src/timetree/inference/forward_pass.rs) commits `max(likely_time, parent_time)` for non-root, non-leaf internal nodes. An inverted mode is therefore stored at exactly the parent's time, permitting a zero-duration edge. Leaves retain their observed dates and are excluded from this projection.

V0 [`ClockTree._ml_t_marginal()`](../../packages/legacy/treetime/treetime/clock_tree.py) commits each independent marginal peak without projection and explicitly notes that marginal reconstruction can produce negative branch lengths. The v1 projection is therefore an intentional-behavior candidate that has not been approved or specified in `kb/decisions/`.

## Inconsistent derived state

The projection changes only the committed node time:

- the node posterior still peaks at the original unconstrained mode;
- confidence extraction computes its highest-posterior-density region from that posterior, then expands the interval to include the projected date;
- positional likelihood evaluates branch distributions at differences between projected node dates;
- each edge's stored `time_length` remains the mode of its independently inferred branch-length distribution.

These values describe different estimators when projection occurs. Equality with the parent also introduces a zero duration without deriving that value from grid spacing, posterior resolution, or date precision.

## Decision required

No implementation ticket is ready until the project chooses the desired contract:

- preserve v0 parity by committing independent marginal modes and handling negative durations explicitly downstream;
- perform constrained inference so posteriors and selected times satisfy topology together; or
- approve a point-estimate projection and define how posterior summaries, confidence intervals, likelihoods, and edge state represent projected nodes.

Validation must include an inverted-mode case and dated-leaf boundaries once the contract is decided.

## Related issues

- [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) - per-pass tail semantics can change posterior peaks but do not jointly constrain marginal modes
