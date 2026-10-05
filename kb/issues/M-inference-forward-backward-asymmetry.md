# GTR transition counting clamps branch lengths and timetree rounds always refit the clock

> [!IMPORTANT]
> **Decision required.** Both instances below reproduce v0 behavior. Changing either one diverges from v0 and needs approval and a `kb/decisions/` entry. The options are:
>
> - Keep v0 parity and close this issue
> - Change the behavior as described under each instance, and accept the divergence from the v0 oracle
>
> Evidence for parity: v0 clamps the branch length used for GTR mutation counting through `_branch_length_to_gtr()` ([packages/legacy/treetime/treetime/treeanc.py#L752-L760](../../packages/legacy/treetime/treetime/treeanc.py#L752-L760), used at [#L1106](../../packages/legacy/treetime/treetime/treeanc.py#L1106)). v0 refits the clock model in every `make_time_tree()` call through `init_date_constraints()` ([packages/legacy/treetime/treetime/clock_tree.py#L418](../../packages/legacy/treetime/treetime/clock_tree.py#L418), [#L374](../../packages/legacy/treetime/treetime/clock_tree.py#L374)), and `TreeTime.run()` calls `make_time_tree()` in every iteration ([packages/legacy/treetime/treetime/treetime.py#L342-L349](../../packages/legacy/treetime/treetime/treetime.py#L342-L349)).

## Instances

### Branch-length floor inside GTR transition counting

The dense and discrete transition counts apply the marginal-pass branch-length floor before they evaluate `expQt` and accumulate dwell times: `fn count_transitions_dense()` calls `effective_branch_length()` ([packages/treetime/src/partition/marginal/shared/data.rs#L32](../../packages/treetime/src/partition/marginal/shared/data.rs#L32), [#L66-L68](../../packages/treetime/src/partition/marginal/shared/data.rs#L66-L68)). The sparse count clamps the same way ([packages/treetime/src/partition/marginal/sparse/count.rs#L29](../../packages/treetime/src/partition/marginal/sparse/count.rs#L29), [#L36](../../packages/treetime/src/partition/marginal/sparse/count.rs#L36)).

GTR parameter estimation therefore sees clamped branch lengths instead of raw values. Zero-length branches contribute dwell time `T_i` at the floor length, which can bias the rate matrix toward short-branch statistics.

Possible change: count transitions with raw branch lengths and handle zero-length branches explicitly.

### Timetree round always refits the clock model

`fn refinement_round()` calls `fn update_clock_model()` unconditionally ([packages/treetime/src/timetree/round.rs#L136](../../packages/treetime/src/timetree/round.rs#L136), [#L385](../../packages/treetime/src/timetree/round.rs#L385)), also when no sequence changed and no polytomy was resolved. The refit is a regression over all dated tips. It costs computation when nothing moved, and small floating-point differences in the regression can perturb the next round and delay convergence detection.

Possible change: skip the refit when the round changed no sequence state, topology, or node time.

## Impact

- GTR inference can be biased by branch-length clamping
- Unconditional clock refits cost computation and can delay convergence detection
