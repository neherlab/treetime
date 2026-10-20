# Config files skip the conflict rules of the command line

clap enforces the conflicts between flags only for values given on the command line. `fn overlay_config` merges the config file into the parsed arguments after clap has run and checks the merged value against the JSON schema only [packages/app-cli/src/cli/config.rs#L19-L44](../../packages/app-cli/src/cli/config.rs#L19-L44). The schema does not carry the clap conflicts, so a config file can set two settings that the command line rejects together.

The app back end already checks these rules for configs: `fn cli_rule_diagnostics()` [packages/app-backend/src/config/cli_rules.rs#L9-L24](../../packages/app-backend/src/config/cli_rules.rs#L9-L24) reads the conflicts and value counts from the clap command and reports each violated rule as a diagnostic. Only the app runs it ([packages/app-backend/src/command.rs#L261](../../packages/app-backend/src/command.rs#L261)). Calling the same check from `overlay_config()` gives the CLI and the app one owner of the rule check; the check can move into `app-commands`, which both surfaces share.

Example: `--keep-root` conflicts with `--reroot` in `clock` ([packages/app-commands/src/commands/clock/args.rs#L140](../../packages/app-commands/src/commands/clock/args.rs#L140)) and in `timetree` ([packages/app-commands/src/commands/timetree/args.rs#L340](../../packages/app-commands/src/commands/timetree/args.rs#L340)).

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
