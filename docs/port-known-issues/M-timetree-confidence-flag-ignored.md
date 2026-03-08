# --confidence flag ignored

`--confidence` alone (with default `--time-marginal=never`) does not produce
`confidence_intervals.tsv`. The flag is declared but never read. The CI output
is gated solely on `time_marginal == OnlyFinal`.

## Root cause

`args.confidence` is declared in
[`args.rs`](../../packages/treetime/src/commands/timetree/args.rs) but never
read anywhere in the timetree pipeline. The CI output is gated on:

```rust
if args.time_marginal == TimeMarginalMode::OnlyFinal {
    // ... write_confidence_intervals ...
}
```

at [`run.rs#L227`](../../packages/treetime/src/commands/timetree/run.rs#L227).

v0 promotes `time_marginal='never'` to `'confidence-only'` (equivalent to
`only-final`) when `--confidence` is set, at
[`wrappers.py#L492`](../../packages/legacy/treetime/treetime/wrappers.py#L492):

```python
time_marginal='confidence-only' if (calc_confidence and time_marginal == 'never') else time_marginal,
```

## v0 covariation requirement

v0 requires `--covariation` for confidence estimation and prints a warning
without it: "Outside of covariation aware mode TreeTime cannot estimate
confidence. Will proceed without confidence estimation." v1 should replicate
this validation rather than silently ignoring both flags.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --confidence \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-confidence data/flu/h3n2/20/aln.fasta.xz
ls tmp/repro-confidence/confidence_intervals.tsv
# Not found (should exist)
```

## Related issues

- [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md) lists
  `--covariation` as another dead flag in the timetree pipeline
