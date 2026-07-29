# Timetree panics on negative numdate during date-string conversion

`year_fraction_to_date()` panics when `year_fraction` is negative (dates before 1 CE). The function uses `f64::fract()`, which preserves the sign, and feeds the result into `std::time::Duration::from_secs_f64()`, which rejects negative values.

Negative numdates arise legitimately when the molecular clock regression extrapolates the root (TMRCA) far into the past. With sparse or low-diversity datasets the clock signal is weak and the extrapolation can overshoot by thousands of years. The value itself is valid math output from the regression; the date-string conversion is what fails.

## Reproduction

```bash
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/lassa/L/20/tree.nwk \
  --dates=data/lassa/L/20/metadata.tsv \
  --output-all=tmp/smoke-tests/timetree/lassa/L/20/basic \
  --output-selection=nwk,nexus,auspice,augur-node-data,clock-model,gtr \
  data/lassa/L/20/aln.fasta.xz
```

Panics with: `cannot convert float seconds to Duration: value is negative`

The offending value: `numdate = -1908.59` (a root extrapolated ~3900 years before the ~2000 CE sampling dates).

## Causal chain

1. Clock regression on 20 Lassa L-segment sequences produces a deep-negative root numdate (`-1908.59`)
2. Augur-node-data output writer calls `year_fraction_to_datestring(numdate)` for every node (`augur_node_data.rs:110`)
3. `year_fraction_to_datestring` delegates to `year_fraction_to_date` (`year_fraction.rs:28-35`)
4. `year_fraction.fract()` returns `-0.59` (sign-preserving for negative inputs)
5. `seconds_in_year as f64 * fraction` produces a negative number of seconds
6. `StdDuration::from_secs_f64(negative)` panics -- `std::time::Duration` is unsigned

## Affected locations

- `packages/treetime-utils/src/datetime/year_fraction.rs:28-35` -- the partial function
- `packages/treetime/src/commands/timetree/output/augur_node_data.rs:110` -- sole crash-triggering call site
- `packages/treetime-utils/src/datetime/parse_date.rs:20` -- other caller, only receives positive input dates (not currently affected)

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
