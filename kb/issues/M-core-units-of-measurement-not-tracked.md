# Units of measurement not tracked or enforced across the codebase

All physical quantities are `f64` with no type-level or documentation-level distinction between units. Functions accept and return bare `f64` values where the caller must know the expected unit from context. This has already caused a concrete bug: `INITIAL_COALESCENT_TC` was set to `0.001` (a subs/site-scale value) when the coalescent model operates in calendar years (PR#851).

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

- `INITIAL_COALESCENT_TC` was `0.001` subs/site when it should have been ~5 years (fixed by PR#851, filed in [N-coalescent-initial-tc-hardcoded-fallback.md](N-coalescent-initial-tc-hardcoded-fallback.md))
- Coalescent Tc coordinate convention (calendar years) not type-enforced ([N-coalescent-time-scale-coordinate-not-type-enforced.md](N-coalescent-time-scale-coordinate-not-type-enforced.md))
- `CalendarTime` newtype exists in coalescent but is not used elsewhere -- node times, date constraints, and confidence intervals pass raw `f64`
- Branch lengths transition from subs/site to years during clock model application with no type change at the boundary
- Optimize command works in subs/site, timetree command works in both subs/site and years depending on pipeline stage -- the transition point is implicit
- `fn ClockModel::clock_deviation(date: f64, div: f64)` accepts interchangeable primitives for calendar date and divergence [`packages/treetime/src/clock/clock_model.rs#L11`](../../packages/treetime/src/clock/clock_model.rs#L11).
- `fn GTR::evolve(t: f64, ...)` calls `t` `branch length (time)`, while lower-level helpers separately distinguish branch-length space after absorbing the model rate [`packages/treetime/src/gtr/gtr.rs#L329`](../../packages/treetime/src/gtr/gtr.rs#L329) [`packages/treetime/src/gtr/gtr.rs#L464`](../../packages/treetime/src/gtr/gtr.rs#L464).
- `TimeLength`, divergence, branch length, rate multiplier, variance, and date APIs all expose raw `f64`, so signatures cannot prevent swapping values with incompatible units.

## Investigation scope

- Audit all command pipelines (`ancestral`, `clock`, `timetree`, `mugration`, `optimize`, `prune`) for implicit unit transitions
- Identify function boundaries where units change without type change
- Evaluate introducing newtype wrappers (`SubsPerSite(f64)`, `CalendarYears(f64)`, `SubsPerSitePerYear(f64)`) at high-risk boundaries
- Assess whether `CalendarTime` from coalescent should be promoted to a shared type
- Document unit conventions per command if type enforcement is deferred

## Related issues

- [N-coalescent-time-scale-coordinate-not-type-enforced.md](N-coalescent-time-scale-coordinate-not-type-enforced.md)
- [N-coalescent-initial-tc-hardcoded-fallback.md](N-coalescent-initial-tc-hardcoded-fallback.md)
