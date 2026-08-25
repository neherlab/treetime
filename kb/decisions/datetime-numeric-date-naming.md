# Numeric date representation and naming

## Decision

v1 uses **numeric date** as the canonical term for a calendar date encoded as a floating-point year: the integer part is the year and the fractional part is the elapsed fraction of that year (for example, `2015.3` is early-to-mid April 2015). The shared scalar type is `DateNumeric` and the shared interval type is `DateRangeNumeric`.

## Context

A single floating-point year value (`YYYY.F`) is the primary time coordinate for molecular-clock and phylodynamic inference. It appears at tree tips (sample dates), on internal nodes (inferred dates), and on coalescent segment boundaries.

The wider field has no single universal term. Two conventions dominate:

- **numeric date** -- BEAST/BEAUti, Nextstrain, augur, and TreeTime v0.
- **decimal date** -- the R ecosystem (`lubridate::decimal_date()`, `ape`).

v0 is consistent: it names the value `numdate`, constructs it with `numeric_date()`, and documents the `YYYY.F` format in [packages/legacy/treetime/treetime/utils.py](../../packages/legacy/treetime/treetime/utils.py).

v1 had drifted into three names for the same value:

- `numdate` -- augur JSON input and output ([packages/treetime/src/commands/timetree/output/augur_node_data.rs](../../packages/treetime/src/commands/timetree/output/augur_node_data.rs)), inherited from the Auspice `num_date` node attribute.
- `CalendarTime` -- the coalescent time coordinate ([packages/treetime/src/coalescent/time_coordinate.rs](../../packages/treetime/src/coalescent/time_coordinate.rs)).
- `year_fraction` -- the datetime utilities ([packages/treetime-utils/src/datetime/year_fraction.rs](../../packages/treetime-utils/src/datetime/year_fraction.rs)).

`year_fraction` is inaccurate: the year fraction is only `F` (the `0.3`), not the whole `2015.3`. `CalendarTime` is ambiguous, because "calendar time" usually denotes an actual calendar date rather than the decimal-year encoding.

## Rationale

- **Ecosystem parity.** v0, augur, and the Auspice schema (`num_date`) already use "numeric date". Matching that term keeps v1 aligned with the tools it reads from and writes to, and avoids a competing vocabulary at the I/O boundary.
- **Single internal vocabulary.** One term (`DateNumeric` / `DateRangeNumeric`) replaces the `numdate` / `CalendarTime` / `year_fraction` split, so the same concept reads the same way across crates.
- **Accuracy over the retired names.** "Numeric date" describes the whole value correctly, unlike "year fraction", and is unambiguous about representation, unlike "calendar time".

## Alternatives considered

- **decimal date** (`DateDecimal` / `DateRangeDecimal`) -- the most accurate literal description and the strongest tooling precedent (`lubridate`, `ape`). Rejected because v1's I/O contracts (augur `num_date`, v0 `numdate`) already fix the ecosystem term as "numeric", and a second term at the boundary would reintroduce the inconsistency this decision removes.
- **Keeping `year_fraction` / `CalendarTime`** -- rejected: `year_fraction` misnames the value and `CalendarTime` is ambiguous.

## Consequences

- Shared and new code names the scalar `DateNumeric` and the interval `DateRangeNumeric`.
- `DateRangeNumeric` is the common interval type for numeric-date spans, including coalescent skyline segments.
- The retired internal names (`year_fraction`, `CalendarTime`) for this value are superseded by the canonical names.
- `numdate` / `num_date` remain only as external field names in the augur/Auspice JSON schema, which is an external contract and is not renamed.
