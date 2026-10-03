# Timetree panics on negative numdate during date-string conversion

`year_fraction_to_date()` panics when `year_fraction` is negative (dates before 1 CE). The function uses `f64::fract()`, which preserves the sign, and feeds the result into `std::time::Duration::from_secs_f64()`, which rejects negative values.

Negative numdates arise legitimately when the molecular clock regression extrapolates the root (TMRCA) far into the past. With sparse or low-diversity datasets the clock signal is weak and the extrapolation can overshoot by thousands of years. The value itself is valid math output from the regression; the date-string conversion is what fails.

## Reproduction

A fixed clock rate far below the substitution rate of the data places the root thousands of years before the samples:

```bash
treetime timetree --tree=data/zika/20/tree.nwk --dates=data/zika/20/metadata.tsv \
  --aln=data/zika/20/aln.fasta.xz --clock-rate=1e-6 --output-all=<dir>
```

The `timetree/zika/20/clock-rate-1e-6` row of `dev/smoke.toml` runs this command. On `data/lassa/L/20` without a fixed rate, the run now stops earlier, at the NaN time distribution of [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md).

Panics with: `cannot convert float seconds to Duration: value is negative`

The offending value is the negative `numdate` of the root.

## Causal chain

1. The clock model places the root before 1 CE, so its `numdate` is negative
2. Augur-node-data output writer calls `year_fraction_to_datestring(numdate)` for every node ([packages/app-output/src/augur_node_data.rs#L97](../../packages/app-output/src/augur_node_data.rs#L97))
3. `year_fraction_to_datestring` delegates to `year_fraction_to_date` ([packages/treetime-utils/src/datetime/year_fraction.rs#L27-L29](../../packages/treetime-utils/src/datetime/year_fraction.rs#L27-L29))
4. `year_fraction.fract()` returns `-0.59` (sign-preserving for negative inputs)
5. `seconds_in_year as f64 * fraction` produces a negative number of seconds
6. `StdDuration::from_secs_f64(negative)` panics -- `std::time::Duration` is unsigned

## Affected locations

- [packages/treetime-utils/src/datetime/year_fraction.rs#L45-L52](../../packages/treetime-utils/src/datetime/year_fraction.rs#L45-L52): `fn year_fraction_to_date`, the partial function. Its lint suppression assumes a within-year second span, but a negative input gives a negative span
- [packages/app-output/src/augur_node_data.rs#L97](../../packages/app-output/src/augur_node_data.rs#L97): the node-data writer, which converts every node date
- [packages/app-commands/src/results/year_date.rs#L18](../../packages/app-commands/src/results/year_date.rs#L18): the run results of the app, which convert dates the same way
- [packages/treetime-utils/src/datetime/parse_date.rs#L28](../../packages/treetime-utils/src/datetime/parse_date.rs#L28): converts input dates, which are positive

## v0 reference behavior

v0 does not crash. `datestring_from_numeric()` (`packages/legacy/treetime/treetime/utils.py:196-214`) catches the exception from `datetime_from_numeric` and falls back to `floor(numdate)` for the year and Python's euclidean `numdate % 1` (always non-negative) for the day-of-year, producing strings like `-1909-05-30`.

## Possible solutions

### A -- fix seam

- **A1 (recommended):** make `year_fraction_to_date` total at the util level. One fix, all current and future callers, mirrors where v0 puts its fallback.
- A2: guard only at the augur call site. Leaves the util partial.

### B -- semantics for negative numdates

- **B1 (v0 parity):** replicate v0's fallback: `floor(numdate)` for year, euclidean remainder for day-of-year. Produces the same strings as v0 (e.g. `-1909-05-30`). Matches the porting default of exact parity.
- **B2 (principled):** total signed arithmetic producing a genuine proleptic Gregorian date. Cleaner than v0's quirky `1900 + frac` reconstruction, but diverges from v0 output -- requires an intentional-change decision.
- **B3 (omit):** emit `None`/omit the `date` field for non-representable numdates. Simplest; diverges from v0 which always emits a string.

The `date` field is cosmetic -- `augur export v2` reads the numeric `numdate`, not the string.

### C -- separate question (not the crash)

Whether `numdate = -1908` for `lassa/L/20` is faithful to v0's clock inference or a v1 clock/extrapolation parity defect. Needs a separate comparison. The crash fix must not mask this.

## Related issues

- [N-timetree-node-data-date-string-fp-boundary.md](N-timetree-node-data-date-string-fp-boundary.md) -- ±1-day rounding in the same `date` field (distinct from this panic)
