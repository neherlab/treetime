# Convert 21 files from loop-over-cases to rstest parameterized

## Problem

Element-by-element loop assertions across 21 test files. Highest density: `test_prop_gtr_site_specific.rs` (39 loops), `generators.rs` (13 loops). Many within proptest bodies where `prop_assert!` is correct but `prop_assert_array_*_eq!` from project utilities could replace the element iteration.

## Files (21)

Listed in source. Exact file list to be determined by searching for the pattern -- the source document states the count as 21 files.

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
