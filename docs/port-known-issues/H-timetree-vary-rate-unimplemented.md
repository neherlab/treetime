# --vary-rate panics with todo!()

The `--vary-rate` flag is accepted by clap but the implementation is a `todo!()`
stub. At runtime, panics with "not yet implemented: calc_rate_susceptibility not
yet implemented".

## Location

- Panic:
  [`run.rs#L224`](../../packages/treetime/src/commands/timetree/run.rs#L224)
- Dead stub:
  [`confidence.rs#L34`](../../packages/treetime/src/commands/timetree/output/confidence.rs#L34)
  `compute_rate_susceptibility()`

## Implementation pieces

Existing infrastructure that a future implementation can use:

- `EdgeTimetree.gamma` (rate multiplier per edge)
- `ClockModel.cov()` (regression covariance)
- `combine_confidence()` (CI combination)

Missing: `compute_rate_susceptibility()` which requires running `run_timetree()`
three times with upper/lower/central rates, and storing per-node date samples for
CI calculation.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --vary-rate \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-vary-rate data/flu/h3n2/20/aln.fasta.xz
# Panics: not yet implemented: calc_rate_susceptibility not yet implemented
```
