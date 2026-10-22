# Timetree convergence metric deficiencies

> [!IMPORTANT]
> **Decision required.** Role of likelihood in `has_converged`: test `log_lh_time`, the log marginal likelihood of the node times; use it as a plateau guard rather than a criterion; or keep the criterion on node times and treat likelihood as reporting only. `log_lh_time` includes the coalescent prior, which is not normalized over node times and changes whenever the coalescent timescale or the lineage counts change, so it is comparable between rounds only while both are unchanged. `log_lh_coal` is the expected complete-data log-likelihood of the current posteriors while the times are inferred under a prior from frozen expected lineage counts (see below), so it is not the maximized objective either

Remaining defects in the timetree convergence tracking: likelihood plays no part in the convergence decision, the reported likelihoods are not comparable across all iterations, and sequence-diff counting undercounts across topology changes.

`fn ConvergenceMetrics::has_converged()` [packages/treetime/src/timetree/convergence/metrics.rs#L26-L31](../../packages/treetime/src/timetree/convergence/metrics.rs#L26-L31) tests node-time movement against `NODE_TIME_TOLERANCE_YEARS` and requires `n_resolved == 0`. It uses `n_diff` only as the fallback when no node is dated on both sides of a round ([kb/decisions/timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md)).

## Details

### Likelihood plays no part in the convergence decision

[packages/treetime/src/timetree/convergence/metrics.rs#L14-L31](../../packages/treetime/src/timetree/convergence/metrics.rs#L14-L31)

`ConvergenceMetrics` carries `log_lh_seq`, `log_lh_time` and `log_lh_coal`, and `has_converged` uses none of them. A tree whose times have settled below tolerance while the likelihood is still moving is declared converged. None of the three is currently a quantity a stop rule could be written against: `log_lh_time` changes with the coalescent timescale and lineage counts (see the decision block above), and `log_lh_coal` is a diagnostic (below).

### The coalescent term is not the objective the times were inferred under

`log_lh_coal` is `fn CoalescentTimescale::expected_log_likelihood()` [packages/treetime/src/timetree/coalescent_timescale.rs#L129](../../packages/treetime/src/timetree/coalescent_timescale.rs#L129): the expected complete-data log-likelihood $-\sum_i I_i/T_{c,i} - \sum_i M_i \ln T_{c,i}$ of the expected lineage counts of the round's posteriors, under the round's $T_c$. The node times are inferred under a prior built from expected lineage counts that stay frozen during the loop and are rebuilt only after a round in which polytomy resolution changed the topology [packages/treetime/src/timetree/round.rs#L54-L58](../../packages/treetime/src/timetree/round.rs#L54-L58) ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)). This is deliberate, because the statistic should describe the current tree, but it makes `log_lh_coal` a diagnostic rather than a term of the maximized objective, and it is not guaranteed to be monotone. Any likelihood-based stop rule has to say which of the two it tests.

### count_sequence_changes underreports on topology changes

`fn count_sequence_changes()` [packages/treetime/src/timetree/convergence/sequence_changes.rs#L9-L13](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L9-L13) sums `fn count_partition_changes()` [packages/treetime/src/timetree/convergence/sequence_changes.rs#L31-L42](../../packages/treetime/src/timetree/convergence/sequence_changes.rs#L31-L42) over the partitions. It compares the internal-node sequences of the two snapshots only for keys present in both. Nodes present in only one snapshot, removed or created by polytomy resolution, are counted as `prev_only` / `curr_only` and logged at debug level, but contribute no diffs. `n_diff` is only the fallback criterion, but the count is still reported in the run log and written to the tracelog. The fix attributes diffs for nodes present in only one snapshot instead of logging and discarding them.

## Impact

- Convergence can be declared while the objective is still moving.
- `log_lh_time` deltas are not usable as a convergence signal across rounds that change the coalescent timescale or lineage counts.
- The tracelog's `n_diff` column undercounts in rounds that changed the topology.

## Validation

Run `data/ebola/20` and `data/mpox/clade-ii/1000`, each with and without `--coalescent opt` and `--keep-polytomies`, and compare the tracelog CSV before and after the change. Round counts must not regress.

## Related

- [N-timetree-convergence-tolerance-vs-branch-grid.md](N-timetree-convergence-tolerance-vs-branch-grid.md)
