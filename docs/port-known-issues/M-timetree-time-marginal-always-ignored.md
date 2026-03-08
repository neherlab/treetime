# --time-marginal=always has no effect

`--time-marginal=always` is accepted by clap but produces the same results as
the default (`never`). v0 runs marginal time reconstruction at every iteration
when `time_marginal='always'`.

## Root cause

`TimeMarginalMode::Always` variant exists in the enum at
[`args.rs#L33`](../../packages/treetime/src/commands/timetree/args.rs#L33) but
is never referenced in the timetree pipeline. Only `TimeMarginalMode::OnlyFinal`
is checked at
[`run.rs#L216`](../../packages/treetime/src/commands/timetree/run.rs#L216) and
[`run.rs#L227`](../../packages/treetime/src/commands/timetree/run.rs#L227).

```
grep -r "TimeMarginalMode::Always" packages/treetime/src/
# No matches
```
