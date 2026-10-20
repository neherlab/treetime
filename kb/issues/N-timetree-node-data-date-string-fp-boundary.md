# Node data `date` string can differ by one day at year-fraction boundaries

The per-node `date` field in timetree node data JSON is the `"YYYY-MM-DD"` string produced from each node's `numdate` by `year_fraction_to_datestring()` [packages/app-output/src/augur_node_data_refine.rs#L129](../../packages/app-output/src/augur_node_data_refine.rs#L129). At year-fraction values just below a day boundary $k/365$ (or $k/366$), the conversion gives the day before the one that v0, and therefore augur, reports from the same `numdate`.

## Cause

v0 writes `node.date = datestring_from_numeric(node.numdate)` [packages/legacy/treetime/treetime/clock_tree.py#L989](../../packages/legacy/treetime/treetime/clock_tree.py#L989). Its `datetime_from_numeric()` adds `1e-10` to the year fraction before it takes the day floor, `int(((numdate % 1) + 1e-10) * days_in_year)`, to avoid this boundary effect [packages/legacy/treetime/treetime/utils.py#L188-L192](../../packages/legacy/treetime/treetime/utils.py#L188-L192).

v1 `fn year_fraction_to_date()` multiplies the fraction by the seconds in the year and adds that duration to January 1, with no such offset [packages/treetime-utils/src/datetime/year_fraction.rs#L40-L47](../../packages/treetime-utils/src/datetime/year_fraction.rs#L40-L47). The `%Y-%m-%d` format then drops the time of day. A fraction that floating-point rounding puts just below $k/365$ therefore lands on day $k - 1$ in v1 and on day $k$ in v0. Example: `numdate = 2019 + 1/365` is `2019.0027397260274` in `f64`; its fraction times the seconds of 2019 is about 86 399.99998 s, so v1 writes `2019-01-01`, while v0 writes `2019-01-02`. For the values $2019 + k/365$, $k = 1, \dots, 364$, the two conversion formulas disagree on 182.

## Impact

The discrepancy is at most one day. `augur export v2` reads `numdate` for the time tree, not the `date` string, so no downstream Auspice output is affected. The mismatch only shows when comparing the `date` field directly against augur's node data JSON.
