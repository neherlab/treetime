# Smoke matrix still runs the retired clock flags, and its snapshots predate the least-squares clock

## Problem

The clock regression became ordinary least squares with an exact root search, and its covariation and optimizer
flags were removed ([clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md),
branch `claude-clock-exact-root`). [dev/smoke.toml](../../dev/smoke.toml) still has rows for the removed flags, which
now fail on argument parsing, and the baseline snapshots hold outputs of the old regression.

Rows to remove (command, variant, arguments):

| command | variant | arguments |
| --- | --- | --- |
| clock | `tip-slack` | `--tip-slack=0.1` |
| clock | `covariation` (flu/h3n2/20, declared to fail) | `--covariation` |
| clock | `covariation` (ebola/20) | `--covariation --sequence-length=19006` |
| clock | `branch-split-brent` | `--branch-split-method=brent` |
| clock | `branch-split-golden-section` | `--branch-split-method=golden-section` |
| clock | `branch-split-grid-21` | `--branch-split-grid-n-points=21` |
| clock | `variance-factor` | `--variance-factor=1` |
| clock | `variance-offset` | `--variance-offset=0.1` |
| timetree | `covariation` | `--covariation` |
| timetree | `tip-slack` | `--tip-slack=0.1` |
| timetree | `reroot-covar` | `--reroot=least-squares --covariation` |
| timetree | `conf-covariation` | `--confidence --covariation --max-iter=3` |

The expected failure of `clock/flu/h3n2/20/covariation` points to `kb/issues/N-clock-regression-all-negative-rate.md`,
which is deleted: the exact root search finds a positive rate where the old search found none.

## Expected output changes

Every case below changes against a baseline from before the change, for the stated reason:

- `clock` and `timetree` with rerooting: the root lies at the exact least-squares optimum of every edge instead of
  an 11-point grid next to the best node; rate, intercept, root position and everything downstream move, by little
  where the old search already found the optimum. Outlier tips no longer enter the fit.
- `optimize --reroot=min-dev`: the root lies at the exact minimum of the root-to-tip variance instead of a Brent
  optimum with tolerance $10^{-6}$; the branch lengths around the root change in the low digits.
- Clock-model JSON: `stats.estimated` has `chisq`, `r_val` and, with a sequence length, `rate_std`; `hessian` and
  `cov` are gone. A timetree with a coalescent prior reports its rate as estimated, not fixed.
- Augur node data of `timetree`: `clock.cov` is gone, `clock.rate_std` is the Gaussian rate standard deviation and
  present whenever the rate is estimated and the sequence length known.
- `--confidence` without `--clock-std-dev`: date intervals include the rate uncertainty (they did only with
  `--covariation` before); with a fixed `--clock-rate` they do not.

## Fix direction

Remove the rows, rebuild the baseline from this branch and record the reasons above with the re-baseline. Check
that no other row relied on the removed flags through a config file (`[config]` cases).

## Validation

`just smoke` (and `just smoke-all`) against the new baseline reports no `changed`, `regressed` or `still-failing`
case.
