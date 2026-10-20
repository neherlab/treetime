# In-loop clock update regresses undated tips on their inferred dates

Inside the refinement loop, v1 re-estimates the clock model from every tip the clock filter did not flag, including tips without a sampling date. `fn update_clock_model()` [packages/treetime/src/timetree/round.rs#L248-L265](../../packages/treetime/src/timetree/round.rs#L248-L265) seeds the regression inputs from `likely_times(graph, constraints, &time_inference)` [packages/treetime/src/timetree/inference/result.rs#L39-L55](../../packages/treetime/src/timetree/inference/result.rs#L39-L55). For a node without a date constraint, `likely_times` returns the peak of the node's posterior time distribution, so an undated tip enters `ClockSet::leaf_contribution` with the date the tree inferred for it.

The fits before the loop read the leaf dates through `fn given_times()` [packages/treetime/src/timetree/inference/result.rs#L28-L37](../../packages/treetime/src/timetree/inference/result.rs#L28-L37) instead: `fn fit_clock_to_dates()` [packages/treetime/src/timetree/optimization/reroot.rs#L75-L83](../../packages/treetime/src/timetree/optimization/reroot.rs#L75-L83), which the initial clock estimate and every pre-loop reroot use, and `fn filter_clock_outliers()` [packages/treetime/src/timetree/pre_loop.rs#L157-L163](../../packages/treetime/src/timetree/pre_loop.rs#L157-L163). There an undated tip has no date and contributes nothing. The rule therefore changes between the pre-loop fits and the in-loop fits.

## Evidence

- **v0**: `ClockTree.setup_TreeRegression` regresses on `tip_value = np.mean(x.raw_date_constraint)` only for terminal nodes with `bad_branch is False` ([packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275)). Undated tips are `bad_branch` and never enter the regression
- **v1 decision text**: [kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md) states that the regression must see the date as given, because regressing on a date the tree inferred feeds the tree's own inference back into the clock it is fitted with. That rule holds for dated tips only; the fallback in `likely_times` bypasses it for undated ones
- **Runtime**: the timetree clock regression table (`--output-clock-csv`) lists the points of the final fit. On `data/zika/20` with the date of one sample removed, the final table after refinement rounds gives that sample a date marked `inferred`, and the sample is not a clock-filter outlier, so it contributes to the fitted rate; with `--max-iter 0` the same sample has no date (test `test_timetree_clock_csv_marks_the_date_of_an_undated_tip` [packages/app-commands/src/commands/timetree/__tests__/test_clock_csv.rs#L99-L104](../../packages/app-commands/src/commands/timetree/__tests__/test_clock_csv.rs#L99-L104), which expects `ClockDateSource::Inferred` after two rounds)

## Impact

- Undated tips pull the fitted clock rate and intercept toward the tree's own time estimates, which diverges from v0 and double counts the tree's information
- The effect grows with the share of undated tips; with every tip dated the fit is unchanged

## Fix direction

[kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md) requires the clock regression to see each leaf date as given. An undated tip has no given date, so under this decision it contributes nothing to the fit, before the loop and inside it. Build the in-loop regression inputs of leaves from the given dates, as the pre-loop fits do, and change the `refined` case of the test above to expect no date for the undated tip.
