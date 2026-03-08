# --relax argument parsing broken

`--relax 1.0 0.5` treats `0.5` as a positional argument (input FASTA path)
instead of the second parameter to `--relax`. The feature is entirely
non-functional.

## Root cause

clap parses `--relax` as taking a single value. The second value `0.5` is
consumed by the positional `[INPUT_FASTAS]...` argument. The `--relax` option
needs to accept exactly 2 values.

## Location

CLI argument definition:
[`args.rs`](../../packages/treetime/src/commands/timetree/args.rs)

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --relax 1.0 0.5 \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --max-iter=3 --outdir=tmp/repro-relax data/flu/h3n2/20/aln.fasta.xz
# Error: When opening file '0.5' - No such file or directory
```
