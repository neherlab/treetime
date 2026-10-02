# In-loop clock update regresses undated tips on their inferred dates

Inside the refinement loop, v1 re-estimates the clock model from every tip the clock filter did not flag, including tips without a sampling date. `fn update_clock_model` in [packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs) seeds the regression inputs from `likely_times(graph, constraints, Some(&inference))` ([packages/treetime/src/timetree/inference/time_inference.rs](../../packages/treetime/src/timetree/inference/result.rs)). For a node without a date constraint, `likely_times` falls back to the peak of the node's posterior time distribution, so an undated tip enters `ClockSet::leaf_contribution` with the date the tree inferred for it.

The steps before the loop (initial estimate, `reroot_tree` in [packages/treetime/src/timetree/optimization/reroot.rs](../../packages/treetime/src/timetree/optimization/reroot.rs)) call `likely_times(.., None)`, where an undated tip has no date and contributes nothing. The rule therefore changes between the pre-loop fits and the in-loop fits.

## Evidence

- **v0**: `ClockTree.setup_TreeRegression` regresses on `tip_value = np.mean(x.raw_date_constraint)` only for terminal nodes with `bad_branch is False` ([packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275)). Undated tips are `bad_branch` and never enter the regression
- **v1 decision text**: [kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md) states that the regression must see the date as given, because regressing on a date the tree inferred feeds the tree's own inference back into the clock it is fitted with. That rule holds for dated tips only; the fallback in `likely_times` bypasses it for undated ones
- **Runtime**: the timetree clock regression table (`--output-clock-csv`) lists the points of the final fit. On `data/zika/20` with the date of one sample removed, the final table after refinement rounds gives that sample a date marked `inferred`, and the sample is not a clock-filter outlier, so it contributes to the fitted rate; with `--max-iter 0` the same sample has no date (test `test_timetree_clock_csv_marks_the_date_of_an_undated_tip` in [`packages/app-commands/src/commands/timetree/__tests__/test_clock_csv.rs`](../../packages/app-commands/src/commands/timetree/__tests__/test_clock_csv.rs))

## Impact

- Undated tips pull the fitted clock rate and intercept toward the tree's own time estimates, which diverges from v0 and double counts the tree's information
- The effect grows with the share of undated tips; with every tip dated the fit is unchanged

## Open question

Whether undated tips should leave the in-loop regression (v0 parity) is a scientific decision for the team.
