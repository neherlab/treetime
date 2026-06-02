# Test build_covariation_clock_params formula

Add unit tests for `timetree::params::build_covariation_clock_params()` which maps covariation options and alignment/sequence length to `ClockParams`.

The function encodes the v0 formula: `variance_factor = 1/seq_len`, `variance_offset_leaf = tip_slack^2 / seq_len^2`, with default `tip_slack = OVER_DISPERSION = 10` (v0 `config.py:8`).

## Test cases

- `covariation=false` returns `None`
- `covariation=true, seq_len=Some(1000), tip_slack=None` -> verify `variance_factor=1e-3`, `variance_offset=0.0`, `variance_offset_leaf=1e-4` (default tip_slack=10)
- `covariation=true, seq_len=Some(500), tip_slack=Some(5.0)` -> verify exact values
- `covariation=true, aln=Some(records)` -> verify seq_len derived from alignment
- `covariation=true, seq_len=None, aln=None` -> error

Expected values from v0 `clock_tree.py:277-285` formula, not from running the SUT.

## Location

`packages/treetime/src/timetree/__tests__/test_params.rs`

## Related issues

Source: `kb/issues/N-test-coverage-gaps.md`
