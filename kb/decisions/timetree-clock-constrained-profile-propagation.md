# Sequence profiles propagate along clock-constrained branch lengths

v1 timetree edges carry two branch lengths. `branch_length` stays the free maximum-likelihood or
input estimate; `clock_branch_length` is the length the inferred node times imply,
$b_e = \mu\,\gamma_e\,(t_{\text{child}} - t_{\text{parent}})$. Marginal reconstruction propagates
profiles along the second. Without this the timetree refinement loop cannot iterate at all.

**Type**: New v1 mechanism reproducing a v0 behavior that v0 obtains without extra work from its choice of
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
[`blended_clock_branch_lengths`](../../packages/treetime/src/timetree/inference/runner.rs#L98),
[`timetree_branch_lengths`](../../packages/treetime/src/timetree/inference/runner.rs#L247-L264), and
the `clock_branch_lengths` map of `struct RoundState` in
[`packages/treetime/src/timetree/round.rs`](../../packages/treetime/src/timetree/round.rs).

## Background

The refinement loop is meant to alternate: reconstruct ancestral states given branch lengths,
re-infer node times given those states, repeat. Without clock-constrained lengths, nothing inside
the loop writes substitution-space branch lengths: `apply_damping`
([`packages/treetime/src/optimize/iteration.rs`](../../packages/treetime/src/optimize/iteration.rs))
writes them only in the branch-length optimization of
[`optimize/pipeline.rs`](../../packages/treetime/src/optimize/pipeline.rs) and
[`optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs) and in the timetree
pre-loop ML step of [`timetree/pre_loop.rs`](../../packages/treetime/src/timetree/pre_loop.rs). The
marginal update was therefore idempotent across rounds and the loop had no way to move. It ran
exactly one iteration on every dataset measured.

The fix is not a free branch-length M-step inside the loop, which would discard the clock. It is
the _constrained_ M-step: profiles move along the lengths the times imply.

The conversion already existed -- `edge_divergence` returns `time_length * rate * gamma` at
[`clock_regression.rs#L474-L486`](../../packages/treetime/src/clock/clock_regression.rs#L474-L486) -- but it
is fed `time_length`, which
[`runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs) sets from
`distribution.likely_time()`, the _unconstrained_ per-edge ML peak, not the inferred
$t_{\text{child}} - t_{\text{parent}}$. So the in-loop clock re-estimation regressed on free ML
lengths and never saw the timetree.

## The two-length split

v0 keeps two lengths deliberately, and v1 now mirrors it:

| role                                                                                                                                                                                         | v0 field          | v1 field                           |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------- | ---------------------------------- |
| ML or input estimate; centres the branch-length grid, feeds the root-to-tip regression (`dist2root` accumulates it, [`treeanc.py#L503`](../../packages/legacy/treetime/treetime/treeanc.py)) | `mutation_length` | `RoundState::branch_lengths`       |
| clock-constrained; feeds profile propagation                                                                                                                                                 | `branch_length`   | `RoundState::clock_branch_lengths` |

v1 had conflated them into one length. Now `blended_clock_branch_lengths()` returns a new map of
clock-constrained lengths, blended with the previous map by the damping factor, and
`timetree_branch_lengths()` builds the per-edge length map for reconstruction: it prefers the
clock-constrained length and falls back to the ML or input length. Only timetree builds this map,
so the ancestral, optimize, and mugration commands are unaffected. It has two production callers:

- `fn RoundState::marginal_update`
  ([`packages/treetime/src/timetree/round.rs#L417-L421`](../../packages/treetime/src/timetree/round.rs#L417-L421)),
  called by `fn refresh_times` (in every `fn refinement_round`) and by `fn final_marginal_round`,
  passes the map to `BranchModel::marginal_update`
  ([`packages/treetime/src/timetree/branch_model.rs#L28-L38`](../../packages/treetime/src/timetree/branch_model.rs#L28-L38))
- `fn assemble_output` passes the map to `fn final_partition_results`
  ([`packages/treetime/src/timetree/pipeline.rs#L150-L158`](../../packages/treetime/src/timetree/pipeline.rs#L150-L158)),
  which runs the final marginal update on the partition set directly

Both reach `PartitionSet::marginal_update`
([`packages/treetime/src/partition/set.rs#L101-L118`](../../packages/treetime/src/partition/set.rs#L101-L118))
and then `MarginalPasses::marginal_update`
([`packages/treetime/src/partition/marginal/shared/update.rs#L48`](../../packages/treetime/src/partition/marginal/shared/update.rs#L48))
through the dense or sparse reconstruction of each partition. Four propagation sites read the map:

- Dense and discrete backward: `fn IndexedBackwardPass::node`
  ([`packages/treetime/src/partition/marginal/shared/pass.rs#L192-L193`](../../packages/treetime/src/partition/marginal/shared/pass.rs#L192-L193))
- Dense and discrete forward: `fn IndexedForwardPass::node`
  ([`packages/treetime/src/partition/marginal/shared/pass.rs#L325`](../../packages/treetime/src/partition/marginal/shared/pass.rs#L325))
- Sparse backward: `fn SparseBackwardPass::node`
  ([`packages/treetime/src/partition/marginal/sparse/backward.rs#L83`](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L83))
- Sparse forward: `fn SparseForwardPass::edge`
  ([`packages/treetime/src/partition/marginal/sparse/forward.rs#L154`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L154))

Transition counting for GTR refinement (`MarginalPasses::count_transitions`) never sees the
clock-constrained lengths: timetree infers its GTR model once, before the loop, from the input
branch lengths (`fn build_partition_set` in
[`packages/treetime/src/partition/set.rs`](../../packages/treetime/src/partition/set.rs), which calls
`fn build_marginal_partition` in
[`packages/treetime/src/partition/create.rs`](../../packages/treetime/src/partition/create.rs) for each file).

## Commit points

`blended_clock_branch_lengths` is called explicitly rather than from inside `run_timetree`, so each
caller chooses its damping and the rate-susceptibility passes can decline to commit at all.

| site                                                                                              | damping                       | why                                                                                                                                                                       |
| ------------------------------------------------------------------------------------------------- | ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `initialize` in `timetree/init.rs`, before the refinement loop                                    | 1.0                           | first commit; nothing to blend with                                                                                                                                       |
| `refine_topology` in `timetree/round.rs`, after polytomy resolution                               | 1.0                           | re-parented edges describe a new pair; the sampled subtree dates every node it creates, so the preliminary times are used directly rather than falling back to ML lengths |
| `refresh_times` in `timetree/round.rs`, after the time inference                                  | `CLOCK_BRANCH_LENGTH_DAMPING` | the loop step                                                                                                                                                             |
| `final_marginal_round` in `timetree/round.rs`, after the final pass of `--confidence` | 1.0                           | reports the final tree                                                                                                                                                    |
| `compute_rate_susceptibility`                                                                     | --                            | does not commit; the three perturbed-rate passes would otherwise leave a blend of the lower-rate and central-rate trees behind                                            |

## Damping

Each pass re-infers every node time at once, which is the classic undamped
alternating-optimization oscillation. With the coalescent lineage counts also held fixed (see
[timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md))
the remaining instability on `data/ebola/20` was a stationary period-2 cycle: the per-round maximum
time change repeated bit-identically at 0.1938358850209 from round 2 onward.

`CLOCK_BRANCH_LENGTH_DAMPING = 0.5` blends half the previously committed length into each commit,
mirroring `apply_damping` in the branch-length optimization loop. The fixed point is unchanged --
$b = (1-f)b + fb$ for any $f$ -- so this changes the path, not the answer.

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
likelihood _falls_, as it must: the clock-constrained lengths are not the ML optimum for the
sequence data, which is the point of a constrained step.

## Related

- [timetree-convergence-on-node-times.md](timetree-convergence-on-node-times.md) -- the loop only
  iterates if the convergence criterion can see the movement this creates
- [M-timetree-polytomy-mutation-counts-read-unconstrained-branch-length.md](../issues/M-timetree-polytomy-mutation-counts-read-unconstrained-branch-length.md)
  -- polytomy mutation counts still read the ML length
