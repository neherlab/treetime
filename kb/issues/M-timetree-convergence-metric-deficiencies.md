# Timetree convergence metric deficiencies

> [!IMPORTANT]
> **Decision required.** Two contracts are open, and no stop rule should be implemented until both are settled:
>
> - Composition of `log_lh_total`: either the total is `None` unless every component is present, or the component set is fixed per run and a missing component is an error. A component dropping out mid-run must not look like a likelihood improvement. `log_lh_coal` goes missing for two reasons: the run has no coalescent (expected), or `collect_coalescent_edges` failed on an inverted edge (a defect, see [M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md))
> - Role of likelihood in `has_converged`: test only the components that are part of the maximized objective; test the total as a plateau guard rather than a criterion; or keep the criterion on node times and treat likelihood as reporting only. The coalescent term is evaluated against live lineage counts while the times are inferred under frozen ones (see below), so the reported total is not the maximized objective

Remaining defects in the timetree convergence tracking system: the reported total log-likelihood is
not comparable across iterations, it plays no part in the convergence decision, and sequence-diff
counting undercounts across topology changes.

The original first defect — the convergence check keying on `n_diff`, which does not measure what
the loop moves — is resolved; see
[timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md).
`has_converged` now tests node-time movement against `NODE_TIME_TOLERANCE_YEARS`, with `n_diff`
retained only as the fallback when no node is dated on both sides of a round.

## Details

### Likelihood plays no part in the convergence decision

[`metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)

`ConvergenceMetrics` carries `log_lh_seq`, `log_lh_pos`, `log_lh_coal` and `log_lh_total`, and
`has_converged` uses none of them. A tree whose times have settled below tolerance while the
likelihood is still moving is declared converged. This cannot be fixed independently of the next
two items: the total is not currently a quantity a stop rule could be written against.

### Total log-likelihood is a sum of different terms across iterations

[`optimizer.rs#L67-L70`](../../packages/treetime/src/timetree/convergence/optimizer.rs#L67-L70)

`log_lh_total` is `[log_lh_seq, log_lh_pos, log_lh_coal].into_iter().flatten().reduce(...)`. When a
component is `Some` in one iteration and `None` in another the number of summed terms changes, so
deltas across iterations compare different objectives. `log_lh_coal` in particular is `None` or
`NaN` whenever `collect_coalescent_edges` fails, which happens silently — see
[M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md).

### The coalescent term is not the objective the times were inferred under

[`likelihood.rs`](../../packages/treetime/src/timetree/convergence/likelihood.rs) evaluates
`compute_coalescent_total_lh` against live lineage counts, while node times are inferred under a
prior built from counts frozen before the loop
([timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)).
This is deliberate — the statistic should describe the tree you actually have — but it means
`log_lh_coal` is a diagnostic rather than a term of the maximized objective, and is not guaranteed
monotone. Any likelihood-based stop rule has to say which of the two it is testing.

### count_sequence_changes underreports on topology changes

[`sequence_changes.rs#L8-L19`](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L8-L19)

Compares per-partition sequence maps by zipping keys present in both snapshots. Nodes present in
only one — removed or created by polytomy resolution — are counted as `prev_only` / `curr_only` and
logged, but contribute no diffs. Lower severity than before, since `n_diff` is now only the
fallback criterion, but the count is still reported and written to the tracelog. The fix
attributes diffs for nodes present in only one snapshot instead of logging and discarding them.

## Impact

- Convergence can be declared while the objective is still moving.
- `log_lh_total` deltas are not usable as a convergence signal without first fixing composition.
- The tracelog's `n_diff` column undercounts in rounds that changed the topology.

## Validation

Run `data/ebola/20` and `data/mpox/clade-ii/1000`, each with and without `--coalescent-opt` and
`--resolve-polytomies`, and compare the tracelog CSV before and after the change. Round counts must
not regress.

## Related

- [N-timetree-convergence-tolerance-vs-branch-grid.md](N-timetree-convergence-tolerance-vs-branch-grid.md)
