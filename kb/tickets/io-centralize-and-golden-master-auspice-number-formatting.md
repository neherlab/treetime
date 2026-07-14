# Centralize and golden-master Auspice number formatting

Use one validated significant-digit formatter and prove parity with the pinned Augur implementation.

## Required changes

- Replace the signed precision parameter with a validated nonzero `u8` fractional-precision type. Reject zero and values outside that type at configuration or deserialization boundaries.
- For each finite nonzero value $n$, define the count of integral digits $d$ and total significant digits $s$ as

  $$
  d = \begin{cases}
  \lfloor \log_{10}(\lfloor |n| \rfloor) \rfloor + 1, & |n| \ge 1, \\
  0, & |n| < 1,
  \end{cases}
  \qquad
  s = d + p,
  $$

  where $p$ is the validated fractional precision. Perform the addition in a wider integer type and return an error when $s > 255$; this is the supported domain of the shared `u8` significant-digit API. Preserve zero unchanged and reject non-finite input, matching the fact that pinned Augur cannot convert non-finite values through its integer-digit calculation and JSON cannot represent them.
- Delete the adapter-local formatting algorithm. Call `float_to_significant_digits(n, s as u8)`, parse its string result back to `f64` with contextual error propagation, and let the Auspice DTO serializer perform the final number-to-JSON conversion. This string-to-number boundary matches Augur's `float(f"{n:.{s}g}")` contract while keeping `treetime-utils` as the sole Rust formatting implementation.
- Capture exact textual expected values by executing pinned Augur only in the project’s trusted reference environment, never from v1.
- Test serialized divergence and date fields, not only helper return values.

## Validation

- Negative values, zero, powers of ten, integer-digit boundaries before and after rounding, tiny magnitudes, ties, non-finite input, and both precision bounds.
- Verify the exact mapping above against pinned Augur revision `d8e38736037ba9474a809f9a5a63bc2b279d2407` for every supported boundary class, including values for which $d + p = 255$ and rejection when $d + p = 256$.
- Debug/release equivalence for the supported domain.
- Whole Auspice JSON golden master and full lint/test suite.

## Related issues

- Source: [kb/issues/N-io-auspice-number-formatting-missing-augur-golden-master.md](../issues/N-io-auspice-number-formatting-missing-augur-golden-master.md)
