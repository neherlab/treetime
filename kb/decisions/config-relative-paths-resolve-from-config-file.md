# Relative paths in config files resolve from the config file's folder

A relative path inside a config file resolves from the folder that contains the config file. This holds for per-command configs (`treetime <command> --config <file>`) and for pipelines (`treetime pipeline --config <file>`).

## Behavior

- **Config file**: a relative input, input template, or output path in the file resolves from the folder of the file
- **Command-line flag**: a path given as a flag resolves from the working directory, also when a config file is given
- **Standard input**: a config read from stdin (`--config -`) has the working directory as its folder
- **Desktop app**: a config file dropped into the app resolves from the dropped file's folder; an example config resolves from its own folder; pasted config text resolves from the working directory of the back end
- **Unchanged values**: `-` (stdin or stdout), absolute paths, and empty strings stay unchanged. `~` is not expanded, and `..` stays in the resolved path
- **Pipelines**: the pipeline file's folder is the base of every step and of the top-level `output_all`. `treetime pipeline --output-all <dir>` sets the top-level output folder from the working directory and wins over the file

v0 has no config files, so this has no v0 counterpart.

## Reason

Example configs must work wherever their folder is unpacked, and a config file must mean the same thing from any working directory. The example configs in `data/` therefore name their inputs relative to their own folder and name no output paths; the user chooses the output folder with `--output-all`.

## Implementation

- `packages/app-commands/src/config/resolve_paths.rs`: the resolver, driven by the `x-path` roles of the command schema
- `packages/app-commands/src/config/load.rs`: per-command config files
- `packages/app-cli/src/cli/pipeline/resolve.rs`: pipelines
- `packages/app-commands/src/check_config.rs`: the config check of the app, with an optional base folder
