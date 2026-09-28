# Residual clock filter computes an exact_date field that is always None

## Problem

`residual_filter()` in `treetime/clock_filter_methods.py` stores `'exact_date': node.raw_date_constraint if type(node) is float else None` for each outlier. `node` is a tree node, never a float, so the value is always `None`. The test was probably meant for `node.raw_date_constraint`.

The outlier table selects only `avg_date`, `tau`, and `residual`, so the field has no effect on the output.

## Fix

Remove the field, or test `node.raw_date_constraint` if a column for exact dates is wanted in the outlier table.
