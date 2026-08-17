# Mass-sized time distributions collapse the forward-pass division to empty

## Symptom and reproduction

With probability-mass grid sizing enabled for timetree distributions
(`rewindow_to_mass` applied at branch-length construction and after the backward/forward composite
steps), the forward pass produces an **empty** time distribution for many internal nodes on some
datasets. Downstream this surfaces as either:

- `Polytomy resolution requires an inferred time for node N, but it has none`, or
- `Coalescent lineage count requires an inferred time for every node, but node (key=N) has none`.

Reproduce by running the timetree pipeline / refinement on a tree that reaches the forward
message-passing division with mass-sized operands (e.g. the polytomy refinement and full-pipeline
integration cases). The golden-master reference runs on other datasets do **not** trigger it, so the
collapse is data-dependent.

## Impact and scope

- Breaks timetree inference (node time assignment) on affected datasets; polytomy resolution and
  coalescent lineage counting then hard-error on the missing times.
- Present whenever mass sizing is active, including construction-only sizing: reverting the
  backward/forward re-windows while keeping only the branch-length construction re-window still
  reproduces the collapse. So the branch-length distributions alone, once mass-sized, are enough to
  trigger it.

## Root cause

The forward refinement computes `parent_except_subtree = distribution_division(parent_time_dist,
msg_to_parent)` to remove a child's own contribution from its parent's posterior. That division is
correct only when the parent posterior and the child message sit on consistent, aligned grids, as
they naturally do when both come from the same message-passing grids.

Mass sizing re-grids a distribution onto its mass-bounded domain (a different extent and spacing) and
does not re-grid `msg_to_parent` to match. The division then evaluates one operand where the other is
at or beyond a hard boundary (`+inf` neg-log = zero probability), so `parent - msg` yields a `-inf`
ordinate (spurious infinite probability) or an all-`+inf` (massless) result. A non-finite peak makes
the subsequent convolution return empty (`convolve.rs` guards a non-finite peak as a massless
operand), and the empty distribution propagates down the tree.

Sizing the exact `[lo, hi]` endpoints (resampling by point count rather than by a `dx` step that can
round a fraction of a cell past a hard bound) removes the `-inf` form at the boundary but not the
massless form: the mass-sized parent's high-probability region can fail to overlap the message grid
where the division is evaluated.

## Workaround

The affected integration tests are marked `#[ignore]` referencing this issue. The mass-sizing feature
otherwise passes its own unit tests (analytic exponential/Laplace domain oracles, re-window mass
conservation, construction mass-fraction and mode-on-hard-bound) and the golden-master reference
runs.

## Fix approach

Make the forward-pass division robust to operands that were re-gridded independently. Options:

- Reconcile the two operands onto a common grid before dividing (resample the parent posterior onto
  the message's support, or both onto a shared fine grid), so the subtraction is always evaluated
  where both carry finite probability.
- Restrict the division result to the finite-mass overlap of the two supports, dropping the
  hard-boundary region where one operand is zero rather than emitting `-inf`/`+inf`.
- Reconsider whether the message-passing posteriors should be stored mass-sized at all, versus mass
  sizing only the final per-node output distributions in a way that never feeds a subsequent
  division (noting that iterative refinement reuses stored distributions to rebuild messages, so a
  naive post-pass does not isolate the arithmetic).

Any fix must keep the golden-master reference runs unchanged and restore the ignored integration
tests to passing.
