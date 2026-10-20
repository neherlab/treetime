# Timetree rounds refit the clock model when nothing changed

`fn refinement_round()` calls `fn update_clock_model()` unconditionally [packages/treetime/src/timetree/round.rs#L68](../../packages/treetime/src/timetree/round.rs#L68) [packages/treetime/src/timetree/round.rs#L255](../../packages/treetime/src/timetree/round.rs#L255), also when no sequence changed and no polytomy was resolved. The refit is a regression over all dated tips. It costs computation when nothing moved, and small floating-point differences in the regression can perturb the next round and delay the convergence test. No cost or delay is measured.

> [!IMPORTANT]
> **Decision required.** The behavior reproduces v0, which refits the clock model in every `make_time_tree()` call through `init_date_constraints()` [packages/legacy/treetime/treetime/clock_tree.py#L418](../../packages/legacy/treetime/treetime/clock_tree.py#L418) [packages/legacy/treetime/treetime/clock_tree.py#L374](../../packages/legacy/treetime/treetime/clock_tree.py#L374), and `TreeTime.run()` calls `make_time_tree()` in every iteration [packages/legacy/treetime/treetime/treetime.py#L342-L349](../../packages/legacy/treetime/treetime/treetime.py#L342-L349). The options are:
>
> - Keep v0 parity and close this issue
> - Skip the refit when the round changed no sequence state, topology, or node time, and record the divergence from v0 in `kb/decisions/`
