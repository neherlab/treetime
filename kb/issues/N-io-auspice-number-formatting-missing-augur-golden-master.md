# Auspice number formatting lacks an independent Augur golden master

> [!IMPORTANT]
> The pinned Augur source is available locally at revision `d8e38736037ba9474a809f9a5a63bc2b279d2407`; an independent golden master remains missing.

## Problem

`pub fn format_number()` [packages/treetime-io/src/tree_ir/auspice.rs#L38](../../packages/treetime-io/src/tree_ir/auspice.rs#L38) reimplements Augur's significant-digit behavior for divergence and numeric dates. Its unit tests use manually written expected values attributed to `augur export_v2.format_number`; they do not capture outputs from a pinned Augur revision.

This leaves boundary behavior unverified for negative values, powers of ten, values crossing an integer-digit boundary after rounding, very small magnitudes, and ties affected by Python and Rust formatting differences. TreeIR round-trip tests cannot serve as an independent oracle because they read output produced by the same adapter.

The public `i32` precision API also performs unchecked `significand + precision`. Extreme values can panic in debug, wrap in release, or request disproportionate formatting allocation. Equivalent significant-digit formatting already exists in `treetime-utils`.

## Potential solutions

- O1. Reuse the shared formatting utility, narrow the public precision domain, and pin it to Augur golden output.
- O2. Keep a format-local implementation with a validated precision type and the same oracle. This preserves adapter ownership but duplicates formatting semantics.

## Recommendation

Use O1. Capture boundary outputs from the pinned Augur implementation, make unsupported precision fallible, and use the shared utility as the sole Rust implementation.

For finite nonzero $n$, let $d = \lfloor \log_{10}(\lfloor |n| \rfloor) \rfloor + 1$ when $|n| \ge 1$ and $d = 0$ otherwise. With validated nonzero fractional precision $p$, Augur requests $s = d + p$ total significant digits. Compute $s$ in a wider integer type, accept it only when $s \le 255$, format with `float_to_significant_digits(n, s as u8)`, and parse that string back to the numeric Auspice DTO field. Preserve zero and reject non-finite input. This defines both the shared-utility mapping and the string-to-number boundary that the golden master must verify.

## Required work

- Capture representative and boundary outputs directly from Augur without reimplementing its formatter in the capture path.
- Compare both numeric values and serialized Auspice JSON, including divergence and date precision.
- Require exact Augur textual formatting over the supported precision domain; numeric equivalence alone does not establish serialized-output parity.
- Cover $d + p = 255$ and reject $d + p = 256$, so the `u8` utility boundary cannot wrap or allocate from an unchecked signed precision.

## Locations

- `pub fn format_number()` [packages/treetime-io/src/tree_ir/auspice.rs#L38](../../packages/treetime-io/src/tree_ir/auspice.rs#L38)
- Unit tests [packages/treetime-io/src/tree_ir/__tests__/test_auspice.rs](../../packages/treetime-io/src/tree_ir/__tests__/test_auspice.rs)

## Related KB items

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)
- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
- [kb/issues/N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
