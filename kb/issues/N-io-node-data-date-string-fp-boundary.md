# Node data `date` string can differ by one day at year-fraction boundaries (IO)

The per-node `date` field in timetree node data JSON is the `"YYYY-MM-DD"` string produced from each node's `numdate` by `year_fraction_to_datestring()` [packages/app-output/src/augur_node_data_refine.rs#L127](../../packages/app-output/src/augur_node_data_refine.rs#L127). At year-fraction values just below a day boundary $k/365$ (or $k/366$), the conversion gives the day before the one that v0, and therefore augur, reports from the same `numdate`.

## Cause

v0 writes `node.date = datestring_from_numeric(node.numdate)` [packages/legacy/treetime/treetime/clock_tree.py#L989](../../packages/legacy/treetime/treetime/clock_tree.py#L989). Its `datetime_from_numeric()` adds `1e-10` to the year fraction before it takes the day floor, `int(((numdate % 1) + 1e-10) * days_in_year)`, to avoid this boundary effect [packages/legacy/treetime/treetime/utils.py#L188-L192](../../packages/legacy/treetime/treetime/utils.py#L188-L192).

v1 `fn year_fraction_to_date()` takes the year with `floor`, multiplies the rest of the value by the seconds in that year, and adds this offset to January 1, with no such offset of the fraction [packages/treetime-utils/src/datetime/year_fraction.rs#L41-L53](../../packages/treetime-utils/src/datetime/year_fraction.rs#L41-L53). The date string then drops the time of day. A fraction that floating-point rounding puts just below $k/365$ therefore lands on day $k - 1$ in v1 and on day $k$ in v0. Example: `numdate = 2019 + 1/365` is `2019.0027397260274` in `f64`; its fraction times the seconds of 2019 is about 86 399.99998 s, so v1 writes `2019-01-01`, while v0 writes `2019-01-02`. For the values $2019 + k/365$, $k = 1, \dots, 364$, the two conversion formulas disagree on 182.

## Impact

The discrepancy is at most one day. `augur export v2` reads `numdate` for the time tree, not the `date` string, so no downstream Auspice output is affected. The mismatch only shows when comparing the `date` field directly against augur's node data JSON.

The same conversion writes the other date strings of v1, which have the same one-day difference at these boundaries:

- The `date` column of the dates table (`--output-dates`): `fn date_fields()` [packages/app-output/src/dates_table.rs#L51](../../packages/app-output/src/dates_table.rs#L51)
- The calendar day of the dates in the app results: `fn YearDate::new()` [packages/app-backend/src/results/year_date.rs#L18](../../packages/app-backend/src/results/year_date.rs#L18)

## Proposed implementation plan

1. Fix the day of the date string at its owner, `pub fn year_fraction_to_datestring()` [packages/treetime-utils/src/datetime/year_fraction.rs#L32](../../packages/treetime-utils/src/datetime/year_fraction.rs#L32): take the day as `floor((fraction + 1e-10) * days_in_year)` as v0 does, and add that many days to January 1. Leave `fn year_fraction_to_date()` unchanged, because `fn year_fraction_days_between()` and other callers need the time of day. All three writers (`augur_node_data_refine.rs`, `dates_table.rs`, `year_date.rs`) call the string function and get the fix with no change
2. Tests in the `__tests__/` of `packages/treetime-utils/src/datetime/`: an `rstest` table with `2019 + 1/365` → `2019-01-02`, `2019 + 364/365` → `2019-12-31`, `2020 + 60/366` → `2020-03-01` (leap year), `2019.0` → `2019-01-01`, and a value just above a boundary (oracle: v0 `datestring_from_numeric()` in `packages/legacy/treetime/treetime/utils.py`, run by hand on the same values). A property test over `k = 1..364` that `2019 + k/365` maps to day `k + 1` (oracle: the calendar)
3. Validate: `just lint-rs`, `just test-rs`, and `just smoke` on a timetree dataset, where only `date` strings at boundaries may change
