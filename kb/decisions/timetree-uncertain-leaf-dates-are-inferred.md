# Uncertain and missing leaf dates are inferred, exact ones are not

The forward pass refines the time distribution of every node whose date is not exact, leaves
included. A leaf given a date range keeps the range as its support and is placed within it by the
message coming down from its parent; a leaf given no date at all is placed by that message alone. A
node given an exact date is left at that date: not refined, and not projected onto its parent's
committed time.

**Type**: v0 parity for uncertain leaves (v0 `ClockTree._ml_t_marginal()` skips only `Delta`-
constrained leaves). New in v1 for undated leaves, which v0 leaves undated.

**v1 location**:
[`refine_distribution_from_parent`](../../packages/treetime/src/timetree/inference/forward_pass.rs)
and the constraint lift in
[`propagate_distributions_backward_slot`](../../packages/treetime/src/timetree/inference/backward_pass.rs).
The input date is stored by
[`load_date_constraints`](../../packages/treetime/src/clock/date_constraints.rs).

## The distinction is exactness, not leafhood

The forward pass used to return early for every leaf, with the comment "do not overwrite leaf
time_distribution (date constraint)". That protects the input, but it also throws away everything
the tree knows about a leaf whose date is uncertain: such a leaf was committed at the midpoint of
its range, whatever the rest of the tree implied, and its confidence interval was the full range.

An uncertain date is not an observation of a time, it is an observation of an interval. The interval
is a factor in the leaf's posterior, exactly as a child message is a factor in an internal node's
posterior, and the two multiply. Only an exact date makes the posterior degenerate, and only then is
there nothing left to infer.

The same rule governs the projection onto the parent time (see
[M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md)):
an inferred node time is clamped to be no earlier than its parent's, an exact date is not. Clamping
an observed date would hide the conflict between that date and the fitted clock that
`commit_clock_branch_lengths` reports as an inverted branch.

## The input date has to be stored separately

The forward pass refines `time_distribution` in place, and the passes run repeatedly: once per
refinement round, plus three more times per rate-susceptibility analysis. The input date is
therefore held in its own field, `NodeTimetree::date_constraint`, written once when the dates are
loaded and never again. Each backward pass multiplies it back into the node's distribution before
building the message to the parent.

Without that, a refined leaf distribution would be sent back up to its parent on the next round as
if it were an independent observation, counting the parent's own message toward the leaf a second
time — an echo that tightens every round while looking like convergence.

`ClockNode::likely_time()`, which the clock regression and the clock filter read for each leaf,
reads the constraint in preference to the time distribution. The regression must see the date as
given: regressing on a date the tree inferred would feed the tree's own inference back into the
clock it is fitted with. Measured on `data/ebola/20`, the fitted clock rate is unchanged to all
printed digits by this change.

## Results

`data/ebola/20` with two dates coarsened to the year (`2014-XX-XX`, `2015-XX-XX`), marginal branch
lengths, `--time-marginal only-final`:

| leaf             | true date | before        | after         | 90% CI after        |
| ---------------- | --------- | ------------- | ------------- | ------------------- |
| `J0037`          | 2014.75   | 2014.50       | 2014.75       | [2014.57, 2015.00]  |
| `0205_C2_PL6086` | 2015.28   | 2015.50       | 2015.14       | [2015.00, 2015.58]  |

Before, both sat at the midpoint of their year with a confidence interval spanning the whole year,
which is what the input said and nothing more. After, they sit where the rest of the tree puts them
within that year. With every date exact, the output is byte-identical to before.

An undated leaf is the limiting case of the same rule: on the run above with `V633`'s date also
removed, it is now dated 2014.55 from its parent and branch-length distribution (true date 2014.61),
where it was previously left undated with a zero-length branch. It stays a bad branch throughout —
it sends no message to its parent and is excluded from the coalescent — so it still contributes
nothing to the fit.

## Related

- [M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md)
  — the projection this rule selects nodes for
- [M-clock-date-interval-collapsed-to-mean.md](../issues/M-clock-date-interval-collapsed-to-mean.md)
  — the clock regression still weights an interval as if it were its midpoint
- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md)
  — the tail policy the forward message carries into the leaf's posterior
