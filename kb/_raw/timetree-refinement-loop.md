# The timetree refinement loop: why it did not iterate, and what makes it converge

Session notes, 2026-08-01. Investigation into streamlining the timetree pipeline
([`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs)), the
measurements that came out of it, and the changes made in response. Measured on `data/ebola/20`
(39 nodes) and `data/mpox/clade-ii/1000` (1917 nodes).

## The starting question

The pipeline is meant to work as an EM-style alternation: rough node times first, then a coalescent
prior estimated from those times, then times re-inferred under that prior, iterating until both
settle. Rerooting and polytomy resolution sit inside the same loop, and a rate-susceptibility
analysis follows it.

## What the pipeline actually did

Three measurements, all taken before any change:

**The loop ran exactly one iteration.** Nothing inside it wrote substitution-space branch lengths
(only the two pre-step calls to `apply_damping` in
[`optimize/iteration.rs`](../../packages/treetime/src/optimize/iteration.rs) do), so
`update_marginal` was idempotent across rounds, so `n_diff == 0` at the end of round 1, so
`has_converged()` fired. Consequences: the `i >= 2` gate on Tc re-estimation was unreachable, **Tc
was estimated exactly once**, and the times ↔ Tc alternation never happened. Both `ebola/20` and
`mpox/1000` confirmed: `Iteration 1: n_diff=0, n_resolved=0 … [converged]`.

**The reported objective silently changed shape.** On mpox, `log_lh_coal` came back `NaN` because a
marginal time inversion made `collect_coalescent_edges` fail (`child older than parent`), yet
`log_lh_total` was still recorded as the sum of the two surviving terms. The one number that could
have driven convergence was not comparable across iterations, and nothing errored.

**Three of eight inference passes recomputed rate-independent work.** In
`compute_branch_length_distribution`, the 300-point per-edge grid and its likelihood evaluation live
in substitution space; only `time_min`/`time_max` use the clock rate. The rate-susceptibility runs
at rate ± σ therefore redid the whole expensive part to produce what is an axis rescaling. Measured:
8 `run_timetree` passes, ~1.0–2.8 s each, ~11 s total on mpox/1000 with
`--coalescent-opt --confidence --covariation`.

A fourth, found later: the output tree's branch lengths are the peaks of the per-edge branch-length
distributions, not the differences of the inferred dates, so they are mutually inconsistent. In an
ebola nexus output, `NODE_0000000` is dated 2014.14 and its child `EM_079497` 2014.26 (Δt = 0.12)
across a branch written as 0.139.

## The missing feedback, and the units it turns on

`update_marginal` propagates profiles along `edge.branch_length()`. Nothing in the loop updated it,
so the loop had no way to move. The fix is not a free branch-length M-step inside the loop; it is to
close the loop through the times, which is the constrained M-step:

```
b_e  =  mu * gamma_e * (t_child - t_parent)
```

v0 gets this for free because it runs the entire timetree in **divergence units**:
`get_time_before_present(numdate) = (today - numdate) * |clock_rate|`
([`utils.py:85`](../../packages/legacy/treetime/treetime/utils.py)), so
`node.branch_length = node.clock_length` ([`clock_tree.py:925`](../../packages/legacy/treetime/treetime/clock_tree.py))
is dimensionally a no-op, and since `use_mutation_length` is always `False`,
`_branch_length_to_gtr` ([`treeanc.py:751-759`](../../packages/legacy/treetime/treetime/treeanc.py))
hands those clock-constrained lengths to the next reconstruction. v1 works in calendar years, so the
conversion has to be explicit.

Note v0 keeps **two** lengths deliberately: `mutation_length` (ML/input) feeds the grid centre and
the root-to-tip regression (`dist2root` accumulates `mutation_length`,
[`treeanc.py:503`](../../packages/legacy/treetime/treetime/treeanc.py)), while `branch_length`
(clock-constrained) feeds profile propagation. v1 had conflated them into one field.

The conversion already existed in v1 — `edge_divergence` returns `time_length * rate * gamma`
([`clock_regression.rs:319-322`](../../packages/treetime/src/clock/clock_regression.rs)) — but it is
fed `time_length`, which is set from `distribution.likely_time()` *before* belief propagation, i.e.
the unconstrained per-edge ML peak, not the inferred `t_child - t_parent`. So the in-loop clock
re-estimation regressed on free ML lengths and never saw the timetree.

## Changes made

1. **`EdgeTimetree.clock_branch_length`**, committed after every forward pass as `mu * gamma * dt`,
   with `HasBranchLength::profile_branch_length()` (default: `branch_length()`) selecting it for
   marginal reconstruction. The default keeps the ancestral and optimize commands untouched.
   `branch_length` stays the ML/input value, mirroring v0's two-field split.
   Negative durations — a leaf dated earlier than its parent, which the forward pass permits because
   it clamps internal children only — are committed as zero and **reported**, not silently floored.
2. **Invalidation on topology change.** A clock branch length describes one parent-child pair, so
   rerooting (splits, merges, inversions) and polytomy resolution (re-parenting) clear it; the edges
   fall back to ML lengths until the next commit.
3. **The post-reroot re-inference is unconditional.** It was gated on `coalescent_tc.is_some()`,
   which meant that without a coalescent nothing recommitted after the reroot cleared the lengths,
   and the first loop round became a no-op (`max_dt = 0.0`, "converged", `log_lh_seq` identical to
   baseline). The pass exists to restore node times, which has nothing to do with the coalescent.
4. **Convergence on node times.** `ConvergenceMetrics` gains `max_time_change` / `rms_time_change`
   (appended to the tracelog CSV); `has_converged` uses them and falls back to `n_diff` only when no
   node times are comparable. `n_diff` measures the sequence block, which is not what the loop
   moves; `n_resolved == 0` is still required.
5. **k(t) decoupling** — see below.
6. **Damping of the write-back** — see below.

## k(t) has two roles, and only one of them may track the times

With the loop live, the branch-length feedback converged cleanly on its own, but the coalescent
destabilized it. Isolating the blocks on ebola/20 showed the destabilizer is **not** Tc estimation:
holding Tc fixed at a hand-picked constant still failed to settle.

The cause is that k(t) plays two roles:

- it sets the merger rate the coalescent **imposes** on node times (prior);
- it is the statistic Tc is **estimated from**, and the likelihood is evaluated against.

`compute_coalescent_model` derived both from the current times, so the prior was self-referential:
each pass moved the times, which moved k(t), which moved the prior the next pass was inferred under.
The loop chased a receding target.

The fix is to decouple them. `run_timetree` now takes the prior as an argument
(`Option<&CoalescentModel>`) instead of a Tc it turns into a model internally, and the caller
decides which k(t) each role gets:

- **Prior**: `CoalescentModel::with_lineage_counts(k0, Tc)`, where `k0` is frozen when the coalescent
  is first introduced and re-frozen **only on topology change** — a reroot or a resolved polytomy
  changes how many lineages exist, not merely when they merge.
- **Statistic**: `optimize_skyline` and `compute_coalescent_total_lh` keep calling
  `CoalescentPrecomputed::from_graph`, i.e. the live k(t). Unchanged.

A Tc present without frozen counts is rejected rather than silently dropping the prior. As a side
effect the rate-susceptibility runs now hold the prior fixed across the three perturbed rates;
previously each run rebuilt k(t) from its own perturbed times, confounding the sensitivity being
measured with the prior's reaction to it.

## Damping

Freezing k0 removed the drift and left a **stationary period-2 cycle** — with fixed Tc, `max_dt`
repeated bit-identically at 0.1938358850209 from round 2 onward. That is the ordinary undamped
alternating-optimization oscillation (each pass re-infers all times at once), so the write-back
blends half of the previously committed length into each commit,
`CLOCK_BRANCH_LENGTH_DAMPING = 0.5`, mirroring `apply_damping` in the optimize loop. The fixed point
is unchanged (`b = (1-f)b + fb` for any f), so this changes the path, not the answer.

## Results

ebola/20, `max_time_change` per round over 10 rounds:

| configuration | live k(t) | frozen k0 | frozen k0 + damping |
| --- | --- | --- | --- |
| fixed rate, no coalescent | 0.030 → 0.0076, monotone | (unchanged) | converges |
| free rate, no coalescent | 0.023 → 0.0 exactly | (unchanged) | converges, 8 rounds |
| fixed rate, fixed Tc = 2.56 | 0.23 → 0.02 → 0.13, `log_lh_seq` drifting −27480.7 → −27491.5 | exact 2-cycle at 0.19384 | 0.0023 ↔ 0.17, one node |
| free rate, `--coalescent-opt` | 0.24–0.39, no decay, `log_lh_coal` drifting 24.9 → 28.0 | 2-cycle at ~0.27, `log_lh_coal` stable ~27.0 | **converges: 0.325 → … → 4e-12, 8 rounds** |

Under live k(t) with the coalescent on, the loop also walked *downhill* on its own reported total
(−27474.9 → −27487.5), trading sequence likelihood for coalescent prior with nothing arbitrating.

**Residual.** The fixed-Tc case still alternates between 0.0023 and 0.171 years, with
`max/rms = 0.171/0.0274 = 6.24 = sqrt(39)` on a 39-node tree: exactly **one** node flipping between
two positions while everything else is still. Local bistability, not global instability. The
configuration that estimates Tc converges.

**Dates move.** Closing the feedback shifts results: mpox root −0.42 y, mean −0.027 y, max 0.47 y;
ebola root +0.097 y.

**Cost.** mpox/1000 without coalescent: 7.3 s (1 round, before) → 37.5 s (10 rounds, after), ~3.4 s
per round. The loop now does the work it was always supposed to do.

## Open items

**A convergence tolerance below the branch grid cannot be met.** mpox reaches a bit-identical
`log_lh_seq` (−303384.78734447353) from round 4 on, while `max_time_change` keeps jittering at
0.016 / 0.032 / 0.048 / 0.080 — exact integer multiples of a ~0.0160 y quantum (≈5.8 days), the
signature of node times snapping between adjacent points of the 300-point branch grid. The current
`NODE_TIME_TOLERANCE_YEARS = 1e-3` is ~16× below what that dataset can resolve, so such runs always
exhaust `--max-iter`. The tolerance should be derived from the grid step, or the criterion should be
the objective with time movement as a guard.

**mpox with `--coalescent-opt` aborts** at round 2: `Coalescent edge has child older than parent`
from `collect_coalescent_edges`. The inversion is pre-existing (it produced the `NaN` above at
baseline) and lives in the statistic role, which the decoupling deliberately leaves alone; it only
became fatal once round ≥ 2 was reachable. Either tolerate it there or fix it at its origin in the
forward pass, which clamps internal children to their parent but leaves observed leaf dates alone.

**Smoothing k0** was not attempted. With the prior frozen, the step function no longer destabilizes
across rounds, but it would still smooth the prior's log-density *within* a pass and is a plausible
cause of the single bistable node.

**Damping factor.** Flat 0.5, versus the exponential schedule the optimize loop uses
(`0.75^(i+1)`, floor 0.01), which would need the round index threaded into the commit.

**Consumers still reading the wrong length.** `apply_relaxed_clock` compares
`act_len = time_length * clock_rate` against `opt_len = branch_length`
([`relaxed_clock.rs:56-57`](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs)),
but both are ML-derived, so gamma has almost no signal to fit; v0 compares `clock_length` against
`mutation_length`. `polytomy.rs:353` has the same issue. Both want `clock_branch_length`.
`count_transitions` still uses ML lengths, which is unobservable today because GTR inference only
runs before the loop.

**Wider restructuring discussed but not done.** A single tracked objective (accumulating the
normalizers that the BP passes and the branch-length grids currently discard, giving a real marginal
likelihood), caching the substitution-space branch likelihood so rate changes are a pure pushforward,
a state + dirty-flag model replacing the prose-documented ordering constraints, and rate uncertainty
by profile likelihood over mu with a mixture posterior instead of the ±σ / probit / quadrature
chain.
