# Negative clock rate causes assertion failure

Clock regression produces a negative rate when temporal signal is weak or the
tree is badly rooted. The code asserts `clock_rate >= 0.0` and panics instead of
handling the error gracefully.

## Location

[`branch_length_likelihood.rs#L38`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L38)

## Affected datasets

- rsv/a/20: always (weak temporal signal)
- dengue/500: always (clock_rate = -2.24e-1)
- flu/h3n2/20 with `--keep-root`: prevents rerooting which fixes the rate
- flu/h3n2/200 with `--keep-root`: same

Datasets that produce valid rates (ebola/20, zika/20) work with `--keep-root`.

## Root cause

When the temporal signal is weak (few informative date constraints relative to
tree depth), or when the tree root is positioned such that tips with later dates
have shorter root-to-tip distances, the regression slope is negative. v0 handles
this by re-rerooting and retrying, or by falling back to a positive rate
estimate. v1 asserts non-negative and panics.

## Repro

```bash
# keep-root on flu (bad root position)
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --keep-root \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-neg-rate data/flu/h3n2/20/aln.fasta.xz
# assertion failed: clock_rate >= 0.0

# rsv (weak temporal signal)
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/rsv/a/20/tree.nwk \
  --dates=data/rsv/a/20/metadata.tsv \
  --name-column=strain --date-column=num_date \
  --outdir=tmp/repro-neg-rate-rsv data/rsv/a/20/aln.fasta.xz
# assertion failed: clock_rate >= 0.0

# dengue (negative rate = -2.24e-1)
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/dengue/500/tree.nwk \
  --dates=data/dengue/500/metadata.tsv \
  --outdir=tmp/repro-neg-rate-dengue data/dengue/500/aln.fasta.xz
# assertion failed: clock_rate >= 0.0
```

## Related issues

- [Negative clock rate silently stored in input-BL mode](M-timetree-negative-rate-input-bl.md)
  bypasses the assertion through a different code path
- [Zero-length branches cause panic](H-timetree-zero-length-panic.md) blocks
  the same datasets through a different mechanism
