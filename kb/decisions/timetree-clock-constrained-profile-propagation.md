# Sequence profiles propagate along clock-constrained branch lengths

v1 timetree edges carry two branch lengths. `branch_length` stays the free maximum-likelihood or
input estimate; `clock_branch_length` is the length the inferred node times imply,
$b_e = \mu\,\gamma_e\,(t_{\text{child}} - t_{\text{parent}})$. Marginal reconstruction propagates
profiles along the second. Without this the timetree refinement loop cannot iterate at all.

**Type**: New v1 mechanism reproducing a v0 behavior that v0 obtains for free from its choice of
units.

**v0 location**: v0 runs the entire timetree in divergence units.
`get_time_before_present(numdate) = (today - numdate) * |clock_rate|` at
[`packages/legacy/treetime/treetime/utils.py#L85`](../../packages/legacy/treetime/treetime/utils.py),
so `node.branch_length = node.clock_length` at
[`packages/legacy/treetime/treetime/clock_tree.py#L925`](../../packages/legacy/treetime/treetime/clock_tree.py)
is dimensionally a no-op, and since `use_mutation_length` is always `False`,
`_branch_length_to_gtr` at
[`packages/legacy/treetime/treetime/treeanc.py#L751-L759`](../../packages/legacy/treetime/treetime/treeanc.py)
hands those clock-constrained lengths to the next reconstruction.

**v1 location**:
[`HasBranchLength::profile_branch_length`](../../packages/treetime-graph/src/edge.rs),
[`EdgeTimetree::clock_branch_length`](../../packages/treetime/src/payload/timetree.rs), and
[`commit_clock_branch_lengths`](../../packages/treetime/src/timetree/inference/runner.rs).

## Background

The refinement loop is meant to alternate: reconstruct ancestral states given branch lengths,
re-infer node times given those states, repeat. Nothing inside the loop wrote substitution-space
branch lengths — only the two pre-step calls to `apply_damping` in
[`optimize/iteration.rs`](../../packages/treetime/src/optimize/iteration.rs) do — so
`update_marginal` was idempotent across rounds and the loop had no way to move. It ran exactly one
iteration on every dataset measured.

The fix is not a free branch-length M-step inside the loop, which would discard the clock. It is
the *constrained* M-step: profiles move along the lengths the times imply.

The conversion already existed — `edge_divergence` returns `time_length * rate * gamma` at
[`clock_regression.rs#L311-L322`](../../packages/treetime/src/clock/clock_regression.rs) — but it
is fed `time_length`, which
[`runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs) sets from
`distribution.likely_time()`, the *unconstrained* per-edge ML peak, not the inferred
$t_{\text{child}} - t_{\text{parent}}$. So the in-loop clock re-estimation regressed on free ML
lengths and never saw the timetree.

## The two-length split

v0 keeps two lengths deliberately, and v1 now mirrors it:

| role | v0 field | v1 field |
| --- | --- | --- |
| ML or input estimate; centres the branch-length grid, feeds the root-to-tip regression (`dist2root` accumulates it, [`treeanc.py#L503`](../../packages/legacy/treetime/treetime/treeanc.py)) | `mutation_length` | `branch_length` |
| clock-constrained; feeds profile propagation | `branch_length` | `clock_branch_length` |

v1 had conflated them into one field. `profile_branch_length()` is a defaulted method on
`HasBranchLength` returning `branch_length()`, overridden only by `EdgeTimetree`, so the
ancestral, optimize, and mugration commands are unaffected.

Six propagation sites read it —
[`marginal_core.rs`](../../packages/treetime/src/partition/marginal_core.rs) (dense backward,
dense forward, indexed backward, indexed forward) and
[`marginal_passes.rs`](../../packages/treetime/src/partition/marginal_passes.rs) (sparse backward,
sparse forward). `count_transitions` deliberately still reads the ML length: GTR inference only
runs before the loop, so the choice is unobservable today and changing it would be an unmeasured
behavior change.

## Commit points

`commit_clock_branch_lengths` is called explicitly rather than from inside `run_timetree`, so each
caller chooses its damping and the rate-susceptibility passes can decline to commit at all.

| site | damping | why |
| --- | --- | --- |
| pipeline, before the optimization loop | 1.0 | first commit; nothing to blend with |
| `Refinement::refine_topology`, after polytomy resolution | 1.0 | re-parented edges describe a new pair; the sampled subtree dates every node it creates, so the preliminary times are used directly rather than falling back to ML lengths |
| `Refinement::run`, after `rebuild_inference` | `CLOCK_BRANCH_LENGTH_DAMPING` | the loop step |
| pipeline, after the final `TimeMarginalMode::OnlyFinal` pass | 1.0 | reports the final tree |
| `compute_rate_susceptibility` | — | does not commit; the three perturbed-rate passes would otherwise leave a blend of the lower-rate and central-rate trees behind |

## Damping

Each pass re-infers every node time at once, which is the classic undamped
alternating-optimization oscillation. With the coalescent lineage counts also held fixed (see
[timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md))
the remaining instability on `data/ebola/20` was a stationary period-2 cycle: the per-round maximum
time change repeated bit-identically at 0.1938358850209 from round 2 onward.

`CLOCK_BRANCH_LENGTH_DAMPING = 0.5` blends half the previously committed length into each commit,
mirroring `apply_damping` in the branch-length optimization loop. The fixed point is unchanged —
$b = (1-f)b + fb$ for any $f$ — so this changes the path, not the answer.

## Backwards branches

A negative duration is committed as zero and reported. The forward pass clamps internal children
to their parent but leaves observed leaf dates as given, so a leaf dated before its parent reaches
the commit; that is a real conflict between the observed date and the fitted clock, not a rounding
artifact, and it is surfaced rather than silently floored. Four such branches occur on
`data/mpox/clade-ii/1000`. The underlying projection contract is
[M-timetree-marginal-node-times-can-violate-topology](../issues/M-timetree-marginal-node-times-can-violate-topology.md).

## Consequences

Closing the feedback changes results. On `data/ebola/20` with no coalescent, switching profile
propagation from ML to clock-constrained lengths moved 12 of 39 node dates by up to 0.0077 y and
changed 25 branch lengths, with `log_lh_seq` going from -27467.34 to -27481.19. The sequence
likelihood *falls*, as it must: the clock-constrained lengths are not the ML optimum for the
sequence data, which is the point of a constrained step.

## Related

- [timetree-convergence-on-node-times.md](timetree-convergence-on-node-times.md) — the loop only
  iterates if the convergence criterion can see the movement this creates
- [M-timetree-consumers-read-unconstrained-branch-lengths.md](../issues/M-timetree-consumers-read-unconstrained-branch-lengths.md)
  — `apply_relaxed_clock` and polytomy mutation counts still read the ML length
