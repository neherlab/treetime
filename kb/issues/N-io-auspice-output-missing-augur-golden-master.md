# Auspice output lacks an independent Augur golden master

No test compares a TreeTime Auspice JSON document, or the numbers in it, with output of a pinned Augur revision. The tests of the Auspice writer check schemas, colorings, and hand-written number examples, and the output-generating tests read output of the same code, so they cannot serve as an independent oracle.

## Number formatting

`pub(crate) fn format_number()` [packages/app-output/src/auspice.rs#L353-L367](../../packages/app-output/src/auspice.rs#L353-L367) reimplements Augur's significant-digit behavior (`export_v2.format_number`) for divergence, numeric dates and their confidence intervals, trait confidences, and trait entropy. The callers pass the fractional precision 6 for divergence and 3 for the other fields [packages/app-output/src/auspice.rs#L304-L338](../../packages/app-output/src/auspice.rs#L304-L338). Its unit tests in [packages/app-output/src/__tests__/test_tree_output.rs](../../packages/app-output/src/__tests__/test_tree_output.rs) use hand-written expected values attributed to Augur; they do not capture output from a pinned Augur revision.

This leaves the boundary behavior unverified for negative values, powers of ten, values that cross an integer-digit boundary after rounding, very small magnitudes, and ties where Python and Rust formatting can differ. Equivalent significant-digit formatting already exists in `treetime-utils`: `float_to_significant_digits()` [packages/treetime-utils/src/fmt/float.rs#L30](../../packages/treetime-utils/src/fmt/float.rs#L30).

## Whole document

The Auspice writer already emits typed nucleotide mutations, genome annotations from the root sequence and the CDS annotations, the `entropy` panel [packages/app-output/src/auspice.rs#L101-L113](../../packages/app-output/src/auspice.rs#L101-L113) [packages/app-output/src/auspice.rs#L160-L184](../../packages/app-output/src/auspice.rs#L160-L184), and `num_date` with its confidence interval for each node [packages/app-output/src/auspice.rs#L191-L226](../../packages/app-output/src/auspice.rs#L191-L226). Whether the complete document, with annotations, panels, mutations, and dates, matches what Augur writes for the same inputs is not checked.

Two Auspice gaps have their own issues: branch support ([M-io-branch-support-dropped-from-outputs.md](M-io-branch-support-dropped-from-outputs.md)) and the `confidence` and `bad_branch` attributes of `ancestral` ([N-ancestral-auspice-missing-v0-attributes-and-default.md](N-ancestral-auspice-missing-v0-attributes-and-default.md)).

## Recommendation

- Use `float_to_significant_digits()` as the only Rust implementation of the formatting. For finite nonzero $n$, let $d = \lfloor \log_{10}(\lfloor |n| \rfloor) \rfloor + 1$ when $|n| \ge 1$ and $d = 0$ otherwise. With the fractional precision $p$, Augur requests $s = d + p$ total significant digits: format with `float_to_significant_digits(n, s)`, parse the string back to `f64`, and let the Auspice serializer write the number. This matches Augur's `float(f"{n:.{s}g}")`
- Capture representative and boundary values, and a whole document for a small dataset, by running a pinned Augur revision in the project's trusted reference environment; never derive expected values from v1 output

## Validation

- Negative values, zero, powers of ten, integer-digit boundaries before and after rounding, tiny magnitudes, ties, and non-finite input match the captured Augur values exactly, as text and as numbers
- A whole-document golden master of the Auspice JSON for a small dataset, with and without genome annotations
- Auspice substitution, insertion, and deletion cases for the verified mutation mappings

## Related KB items

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
