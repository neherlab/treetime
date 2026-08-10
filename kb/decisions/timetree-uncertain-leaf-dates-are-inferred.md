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

## Why a posterior can come out empty at all

It should not be able to. Every constraint is already in the backward pass, and the forward pass
only picks the likely time within them, so write the node's posterior as

```
P(t) = cavity(t)  x  subtree(t),      cavity = (parent posterior / message this node sent up) * branch
```

The parent's posterior is non-empty, and it is itself the product of that cavity and this node's
message. So there is a parent time $t_p$ where both factors are positive, and the message being
positive at $t_p$ means there is a $t$ in this node's own support that the branch can reach from it.
At that $t$ the product above is positive. An empty posterior is therefore never a statement about
the constraints; it is always an artifact of how the distributions are represented.

Three ways the representation produces one, in falling order of how often they are hit:

1. **Underflow.** The distributions are gridded plain probabilities. A posterior that is merely
   astronomically small evaluates to `0.0`, and `normalize()` maps a maximum of zero to `Empty`.
   This is [M-timetree-backward-pass-plain-space-underflow.md](../issues/M-timetree-backward-pass-plain-space-underflow.md)
   seen from the forward side.
2. **Truncation.** A grid ends where the values stopped mattering, and a `Hard` tail past its end
   states the value there is zero rather than small. Two supports that overlap in exact arithmetic
   can then be declared disjoint.
3. **Collapse to a point.** When two supports meet at a single grid point the product is a `Point`.
   A delta annihilates every range that does not contain it, so a whole subtree of dates can be
   contradicted at once.

The two nodes contradicted on `data/mpox/clade-ii/2000` are case 1, and the underflow is telling the
truth. Both are clock-filter outliers whose divergence implies a date around 2110-2120 against a
stamped 2022; the message reaching them spans `[2003.9, 2477.1]` on 313 points, so at their stamped
date the density is far below what an `f64` holds. Keeping the given date is the right answer for
such a tip, and the clock filter reports the same sequences separately.

## A date the tree contradicts is kept, not refined away

The product of the message from the parent and a date range is empty when the two are disjoint, and
a node refined onto an empty distribution is left undated. That is strictly worse than the date the
input carries, so the refinement is abandoned for such a node and the given date stands. Nodes are
counted and the count reported once per pass; `--verbosity=debug` names them and prints what the
tree implied for each. A node whose distribution comes from its subtree alone is still left undated
by an empty product, which is the older contract and a different situation: there is no input date
there to fall back on.

The disagreement is real often enough to matter. On `data/mpox/clade-ii/2000`, where 629 of 1989
tips carry a month-precision date, two tips are contradicted this way. It is also the failure mode
of a run whose distributions have collapsed onto single points, where every date range in a subtree
is contradicted at once; falling back keeps such a run producing the dates it was given rather than
an increasingly undated tree.

## Cost, and the window that pays for it

A node is refined by convolving the message from its parent with the branch-length distribution,
which costs the product of the two grid sizes. Posteriors deep in a large tree carry tens of
thousands of grid points spanning centuries -- 53 000 points over 1733-2019 on
`data/mpox/clade-ii/2000` -- so each refinement is milliseconds, and refining the leaves as well as
the internal nodes made the pass markedly slower.

The product that follows the convolution is zero outside the node's own support, so only the parent
times from which the branch can reach that support contribute anything. `restrict_to_reachable` cuts
the parent to that window, snapped outwards onto the parent's own grid points, before convolving. A
date range spans weeks, so the window is a few hundred points instead of tens of thousands. Internal
nodes benefit too, since a subtree posterior is usually far narrower than the parent's full support.

`data/mpox/clade-ii/2000`, `--max-iter 3`, three alternating runs each:

| build                                          | wall clock     | undated tips |
| ---------------------------------------------- | -------------- | ------------ |
| before uncertain leaves were refined           | 65, 66, 68 s   | 0            |
| refining them, without the window or fallback  | 96, 98, 99 s   | 8            |
| refining them, with both                       | 31, 32, 38 s   | 0            |

Exactly dated tips come out identical with and without the window. Uncertain tips move by at most
0.009 y and internal nodes by at most 0.05 y: the window changes which grid the same posterior is
sampled on, and `likely_time` is an argmax over that grid. Both are the size of the quantization the
branch grid already imposes on node times, and below the 1e-2 y the loop calls converged
([timetree-convergence-on-node-times.md](timetree-convergence-on-node-times.md)).

## Related

- [M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md)
  — the projection this rule selects nodes for
- [M-clock-date-interval-collapsed-to-mean.md](../issues/M-clock-date-interval-collapsed-to-mean.md)
  — the clock regression still weights an interval as if it were its midpoint
- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md)
  — the tail policy the forward message carries into the leaf's posterior
- [M-timetree-branch-grid-uniform-resolution.md](../issues/M-timetree-branch-grid-uniform-resolution.md)
  — the grid sizes the reachable window works around rather than fixes
