# --metadata not required, misleading error when omitted

`timetree` declares `--metadata` as optional ([packages/app-commands/src/commands/timetree/args.rs#L142-L149](../../packages/app-commands/src/commands/timetree/args.rs#L142-L149)), and its argument conversion checks only `--tree` ([packages/app-commands/src/commands/timetree/args.rs#L64-L67](../../packages/app-commands/src/commands/timetree/args.rs#L64-L67)). The sibling `clock` command reports a missing `--metadata` with the required-arguments error ([packages/app-commands/src/commands/clock/args.rs#L42-L49](../../packages/app-commands/src/commands/clock/args.rs#L42-L49)).

Without dates, `timetree` runs up to the initial clock fit ([packages/treetime/src/timetree/pipeline.rs#L296-L297](../../packages/treetime/src/timetree/pipeline.rs#L296-L297)) and fails there with a message that does not mention the missing dates:

```bash
./dev/docker/run just r treetime timetree \
  --tree=data/flu/h3n2/20/tree.nwk \
  --alignment=data/flu/h3n2/20/aln.fasta.xz \
  --output-all=tmp/no-dates
```

```
Error:
   0: Failed to infer clock model
   1: Clock rate is negative for all root positions. The data may lack temporal signal. Please specify --clock-rate explicitly.

Location: packages/treetime/src/clock/find_best_root/find_best_root.rs:69
```

With `--keep-root`, the same run fails with "No variation in sampling dates! Please specify your clock rate explicitly." Both messages advise `--clock-rate`, which does not help: with `--clock-rate=0.003`, the run fails at the outlier filter with "Clock filtering requires at least one dated leaf".

The flag should be required for the timetree command, or the error message should suggest providing `--metadata`.
