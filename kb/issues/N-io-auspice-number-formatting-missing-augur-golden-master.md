# Auspice number formatting lacks an independent Augur golden master

> [!IMPORTANT]
> **More research is required.** Current examples support the intended formatting rule, but the claimed Augur parity has not been verified against an independent Augur execution.

## Problem

`pub fn format_number()` [packages/treetime-io/src/tree_ir/auspice.rs#L38](../../packages/treetime-io/src/tree_ir/auspice.rs#L38) reimplements Augur's significant-digit behavior for divergence and numeric dates. Its unit tests use manually written expected values attributed to `augur export_v2.format_number`; they do not capture outputs from a pinned Augur revision.

This leaves boundary behavior unverified for negative values, powers of ten, values crossing an integer-digit boundary after rounding, very small magnitudes, and ties affected by Python and Rust formatting differences. TreeIR round-trip tests cannot serve as an independent oracle because they read output produced by the same adapter.

## Research required

- Identify and pin the Augur revision that defines the compatibility contract.
- Capture representative and boundary outputs directly from Augur without reimplementing its formatter in the capture path.
- Compare both numeric values and serialized Auspice JSON, including divergence and date precision.
- Decide whether exact Augur formatting is required or whether a documented numeric-equivalence contract is sufficient.

## Locations

- `pub fn format_number()` [packages/treetime-io/src/tree_ir/auspice.rs#L38](../../packages/treetime-io/src/tree_ir/auspice.rs#L38)
- Unit tests [packages/treetime-io/src/tree_ir/__tests__/test_auspice.rs](../../packages/treetime-io/src/tree_ir/__tests__/test_auspice.rs)

## Related KB items

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
- [kb/issues/N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
