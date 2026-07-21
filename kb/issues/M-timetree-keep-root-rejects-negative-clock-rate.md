# Timetree --keep-root rejects negative clock rate

When `--keep-root` is used, rerooting is disabled. If the input tree's root yields a negative clock rate from regression, v1 errors with "Estimated clock rate is non-positive" at [`packages/treetime/src/clock/clock_model.rs#L137`](../../packages/treetime/src/clock/clock_model.rs#L137). v0 proceeds with the negative rate and produces output.

## Reproduction

### v1 (fails)

```bash
# flu/h3n2/20: errors with rate -3.755498e-3
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta.xz \
  --output-all=tmp/timetree/flu-h3n2-20-keep-root \
  --output-selection=nwk,nexus,auspice,augur-node-data,clock-model,gtr \
  --keep-root

# flu/h3n2/200: same error
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/flu/h3n2/200/tree.nwk \
  --dates=data/flu/h3n2/200/metadata.tsv \
  --aln=data/flu/h3n2/200/aln.fasta.xz \
  --output-all=tmp/timetree/flu-h3n2-200-keep-root \
  --output-selection=nwk,nexus,auspice,augur-node-data,clock-model,gtr \
  --keep-root
```

Error: `Estimated clock rate is non-positive (-3.755498e-3)` at `packages/treetime/src/clock/clock_model.rs:137`, called from `ClockRerootResult::into_clock_model()` at `packages/treetime/src/clock/clock_regression.rs:58`, called from `timetree::pipeline::run` at `packages/treetime/src/timetree/pipeline.rs:154`.

### v0 (succeeds)

```bash
# Decompress first -- v0 does not support compressed FASTA
xz -dk data/flu/h3n2/20/aln.fasta.xz

# flu/h3n2/20: succeeds with rate -3.298e-03, r^2=0.74
./dev/docker/python treetime \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta \
  --keep-root \
  --outdir=tmp/v0/flu-h3n2-20-keep-root

# flu/h3n2/200: succeeds with rate -2.085e-03, r^2=0.75
xz -dk data/flu/h3n2/200/aln.fasta.xz
./dev/docker/python treetime \
  --tree=data/flu/h3n2/200/tree.nwk \
  --dates=data/flu/h3n2/200/metadata.tsv \
  --aln=data/flu/h3n2/200/aln.fasta \
  --keep-root \
  --outdir=tmp/v0/flu-h3n2-200-keep-root
```

v0 produces: timetree (nexus, pdf), molecular clock file, ancestral sequences, auspice JSON, divergence tree, dates TSV. The negative rate is recorded in `molecular_clock.txt`.

## v1 behavior

`ClockModel::from_regression()` at `clock_model.rs:136` checks `regression.clock_rate <= 0.0` and returns an error. This check fires during timetree pipeline at `pipeline.rs:154` before any inference begins.

Call chain: `run_timetree_estimation` -> `pipeline::run` -> `ClockRerootResult::into_clock_model` -> `ClockModel::from_regression` -> error.

## Affected datasets

| Dataset      | v1 estimated rate | v0 estimated rate |
| ------------ | ----------------- | ----------------- |
| flu/h3n2/20  | -3.755498e-3      | -3.298e-03        |
| flu/h3n2/200 | (same error)      | -2.085e-03        |

Both are in the smoke test "Known failures" section (`dev/run-smoke-tests` lines 319-320).

## Analysis

The negative rate check is reasonable for the default rerooting path -- if all root positions yield negative rates, the data lacks temporal signal. But with `--keep-root`, the user explicitly opted to keep the input root. Rejecting a negative rate contradicts that intent.

v0's approach (proceed with negative rate) preserves user agency. The rate is reported in output, so the user can evaluate whether results are meaningful.

## Fix direction

When `--keep-root` is active, either:

- Skip the positive-rate check and proceed (v0 parity)
- Warn about negative rate and proceed
- Allow `--clock-rate` to override as the error message suggests, but also allow `--keep-root` to imply acceptance of the regression result

## Related issues

- [N-clock-regression-all-negative-rate.md](N-clock-regression-all-negative-rate.md) -- clock command also rejects all-negative rates, but that path reroots first
