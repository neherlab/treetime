# Timetree rounds refit the clock model when nothing changed

`fn refinement_round()` calls `fn update_clock_model()` unconditionally [packages/treetime/src/timetree/round.rs#L65](../../packages/treetime/src/timetree/round.rs#L65) [packages/treetime/src/timetree/round.rs#L258](../../packages/treetime/src/timetree/round.rs#L258), also when no sequence changed and no polytomy was resolved. The refit is a regression over all dated tips. It costs computation when nothing moved, and small floating-point differences in the regression can perturb the next round and delay the convergence test. No cost or delay is measured.

> [!IMPORTANT]
> **Decision required.** The behavior reproduces v0, which refits the clock model in every `make_time_tree()` call through `init_date_constraints()` [packages/legacy/treetime/treetime/clock_tree.py#L418](../../packages/legacy/treetime/treetime/clock_tree.py#L418) [packages/legacy/treetime/treetime/clock_tree.py#L374](../../packages/legacy/treetime/treetime/clock_tree.py#L374), and `TreeTime.run()` calls `make_time_tree()` in every iteration [packages/legacy/treetime/treetime/treetime.py#L342-L349](../../packages/legacy/treetime/treetime/treetime.py#L342-L349). The options are:
>
> - Keep v0 parity and close this issue
> - Skip the refit when the round changed no sequence state, topology, or node time, and record the divergence from v0 in `kb/decisions/`

## Proposed implementation plan

The cited locations are current: `fn refinement_round()` [packages/treetime/src/timetree/round.rs#L37](../../packages/treetime/src/timetree/round.rs#L37) builds `RoundOutcome` and then calls `fn update_clock_model()` at line 65 without a condition.

- **Keep v0 parity**: delete this issue; no code change
- **Skip the refit when nothing changed**:
  1. Measure first: count the rounds in which `RoundOutcome` reports zero `sequence_changes`, an unchanged `topology`, and a `time_change` of zero on the `data/` datasets, and the time of `fn update_clock_model()` in those rounds. Stop and keep parity when the cost is negligible
  2. Add a method on `RoundOutcome` that says whether the round changed anything, and call `fn update_clock_model()` in `fn refinement_round()` only when it returns true
  3. Test in `packages/treetime/src/timetree/__tests__/`: a round on a converged state returns the clock model unchanged, compared as a whole value with the clock model of the input state
  4. With approval, record the divergence from v0 in `kb/decisions/`, and validate with `just smoke` that node dates stay within `1e-6` of the parity build on the smoke datasets
