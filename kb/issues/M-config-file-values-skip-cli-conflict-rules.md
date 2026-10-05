# Config files skip the conflict rules of the command line

clap enforces the conflicts between flags only for values given on the command line. `fn overlay_config` merges the config file into the parsed arguments after clap has run and checks the merged value against the JSON schema only [packages/app-cli/src/cli/config.rs#L17-L42](../../packages/app-cli/src/cli/config.rs#L17-L42). The schema does not carry the clap conflicts, so a config file can set two settings that the command line rejects together.

Example: `--keep-root` conflicts with `--reroot` in `clock` ([packages/app-commands/src/commands/clock/args.rs#L221](../../packages/app-commands/src/commands/clock/args.rs#L221)) and in `timetree` ([packages/app-commands/src/commands/timetree/args.rs#L439](../../packages/app-commands/src/commands/timetree/args.rs#L439)).

```
treetime clock --tree=data/zika/20/tree.nwk --dates=data/zika/20/metadata.tsv --keep-root --reroot=oldest --output-all=<dir>
error: the argument '--keep-root' cannot be used with '--reroot <REROOT>'
```

The same two settings in a config file run to completion with exit code 0:

```yaml
tree: ../zika/20/tree.nwk
metadata: ../zika/20/metadata.tsv
keep_root: true
reroot: oldest
```

`treetime clock --config=data/smoke/clock-keep-root-and-reroot.yaml --output-all=<dir>` runs it.

## Smoke coverage

The `config-keep-root-and-reroot` CLI row of `dev/smoke.toml` runs `data/smoke/clock-keep-root-and-reroot.yaml` and declares the failure the command line gives. The row reports an unexpected pass until config files apply the same conflicts.
