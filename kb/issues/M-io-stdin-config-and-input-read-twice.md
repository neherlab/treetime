# Config file and an input can both read standard input

A command accepts `--config -` together with `-` for an input, such as `--tree -`. The config file then consumes standard input, and the input reads an empty standard input, so the run fails with an error that names the wrong cause.

## Current state

- `fn check_run_paths()` rejects `-` given for two inputs of one run before any input is read [packages/app-commands/src/commands/shared/run_paths.rs](../../packages/app-commands/src/commands/shared/run_paths.rs). It reads the inputs from the `x-path` roles of the command schema, and `--config` is a CLI-only setting outside the schema, so the check does not see it
- `fn overlay_config()` reads `--config -` from standard input to its end before the run starts [packages/app-cli/src/cli/config.rs](../../packages/app-cli/src/cli/config.rs)

## Reproduction

`treetime prune --config - --tree - --output-tree-nwk out.nwk < config.yaml` reads the config file, then stops with a Newick parse error of the empty standard input.

## Fix direction

In `fn overlay_config()`, when `--config` is `-`, reject a merged config whose inputs contain `-`, with the message of `fn check_run_paths()` that names `--config` and the input.

## Related

- [kb/features/io.md](../features/io.md) -- standard input entry
