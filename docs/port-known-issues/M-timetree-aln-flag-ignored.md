# --aln flag silently ignored

`--aln=file.fasta` is accepted by clap (hidden flag with
`visible_alias("aln")`) but the file is never loaded. The pipeline then fails
with "Alignment required when branch_length_mode is not 'input'" because
`args.input_fastas` is empty.

Users migrating from v0 (which uses `--aln`) will hit this.

## Root cause

`args.aln` is declared at
[`args.rs#L70`](../../packages/treetime/src/commands/timetree/args.rs#L70) but
never read. Alignment loading at
[`initialization.rs#L47`](../../packages/treetime/src/commands/timetree/initialization.rs#L47)
only checks `args.input_fastas`.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta.xz \
  --outdir=tmp/repro-aln
# Error: Alignment required when branch_length_mode is not 'input'
```

## Related issues

- [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md) lists other
  accepted-but-ignored flags
