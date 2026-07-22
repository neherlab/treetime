# Node data `date` string can differ by one day at year-fraction boundaries

The per-node `date` field in timetree node data JSON is the `"YYYY-MM-DD"` string produced from each node's `numdate` by `year_fraction_to_datestring()`. At year-fraction values that fall on a floating-point day boundary, the conversion can round to a day adjacent to the one augur reports from the same `numdate`.

The discrepancy is at most one day. `augur export v2` reads `numdate` for the time tree, not the `date` string, so no downstream auspice output is affected. The mismatch only shows when comparing the `date` field directly against augur's node data JSON.
