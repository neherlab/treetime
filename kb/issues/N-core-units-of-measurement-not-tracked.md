# Units of measurement not tracked or enforced across the codebase

All physical quantities are `f64` with no type-level or documentation-level distinction between units. Functions accept and return bare `f64` values where the caller must know the expected unit from context.

## Known quantities and their units

| Quantity                    | Units                 | Where                                     |
| --------------------------- | --------------------- | ----------------------------------------- |
| Branch length (raw)         | subs/site             | tree input, GTR, optimization             |
| Clock rate                  | subs/site/year        | clock regression                          |
| Node time                   | decimal calendar year | time distributions, forward/backward pass |
| Tc (coalescent timescale)   | years                 | coalescent model                          |
| Branch length (time-scaled) | years                 | after clock model application             |
| Root-to-tip distance        | subs/site             | clock regression input                    |
| Date constraints            | decimal calendar year | metadata input                            |
| Confidence intervals        | years                 | timetree output                           |

## Known problems

- Coalescent Tc coordinate convention (calendar years) not type-enforced ([N-coalescent-time-scale-coordinate-not-type-enforced.md](N-coalescent-time-scale-coordinate-not-type-enforced.md))
- Node times, date constraints, and confidence intervals pass raw `f64` calendar dates
- Branch lengths transition from subs/site to years during clock model application with no type change at the boundary
- Optimize command works in subs/site, timetree command works in both subs/site and years depending on pipeline stage -- the transition point is implicit
- `fn ClockLine::clock_deviation(date: f64, div: f64)` accepts interchangeable primitives for calendar date and divergence [`packages/treetime/src/clock/clock_model.rs#L15`](../../packages/treetime/src/clock/clock_model.rs#L15)
- `fn GTR::evolve(profile, t: f64, return_log)` takes the branch length `t` as a bare `f64` [`packages/treetime/src/gtr/gtr.rs#L128`](../../packages/treetime/src/gtr/gtr.rs#L128)
- Divergence, branch length, rate multiplier, variance, and date APIs all expose raw `f64`, so signatures cannot prevent swapping values with incompatible units

## Investigation scope

- Audit all command pipelines (`ancestral`, `clock`, `timetree`, `mugration`, `optimize`, `prune`) for implicit unit transitions
- Identify function boundaries where units change without type change
- Evaluate introducing newtype wrappers (`SubsPerSite(f64)`, `CalendarYears(f64)`, `SubsPerSitePerYear(f64)`) at high-risk boundaries
- Document unit conventions per command if type enforcement is deferred

## Related issues

- [N-coalescent-time-scale-coordinate-not-type-enforced.md](N-coalescent-time-scale-coordinate-not-type-enforced.md)

## Proposed implementation plan

The cited signatures are current (`clock_deviation()` is a trait method with two implementations at `clock_model.rs` L185 and L208; `GTR::evolve()` is `pub(crate)`), and no unit newtype exists in the workspace.

1. Audit: list every function in `packages/treetime/src/{clock,timetree,coalescent,optimize,gtr}` and `packages/app-commands/src/commands/` that takes or returns two or more `f64` values of different units, and every point where branch lengths change from subs/site to years. Record the list in this issue
2. Add newtypes with `derive_more` (`From`, `Into`, arithmetic derives where the unit allows it) in `packages/treetime-primitives`: `SubsPerSite`, `NumericDate` (the glossary term for calendar dates as `YYYY.F`), `Years` (durations), and `SubsPerSitePerYear`. Give them `serde(transparent)` so no output format changes
3. Apply them first at the high-risk boundaries the audit finds, starting with `ClockLine::clock_deviation(date: NumericDate, div: SubsPerSite)` and the subs/site-to-years conversion of the clock model, then the Tc coordinate of [N-coalescent-time-scale-coordinate-not-type-enforced.md](N-coalescent-time-scale-coordinate-not-type-enforced.md). Keep `ndarray` hot loops on `f64`, and convert at function entry
4. Tests: none for the newtypes themselves (the type checker covers them); existing tests must pass unchanged, which proves the change keeps behavior
5. Validate: `just lint-rs`, `just test-rs`, `just smoke` (no output may change)
