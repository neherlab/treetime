# Config files spell output selections differently from the command line

The per-command output selection enums (`AncestralOutputSelection`, `TimetreeOutputSelection`, `ClockOutputSelection`, `MugrationOutputSelection`, `OptimizeOutputSelection`, `PruneOutputSelection`) derive `Serialize` and `Deserialize` without `#[serde(rename_all = "kebab-case")]`. A config file, a pipeline step, and a server or desktop request must therefore spell them as Rust variant names, while the command line takes the kebab-case names of clap's `ValueEnum`:

- Command line: `--output-selection nwk,mat-pb,augur-node-data`
- Config file: `output_selection: [Nwk, MatPb, AugurNodeData]`; `nwk` is rejected as "not a valid value"

The generated schemas `packages/schemas/input-config-<command>.schema.json` publish the PascalCase values, and so do the `<Command>Config` components of `packages/app-contracts/openapi.json` and the generated zod schemas.

## Evidence

- `macro_rules! per_command_output_selection` [packages/app-commands/src/commands/shared/output_args.rs#L12-L43](../../packages/app-commands/src/commands/shared/output_args.rs#L12-L43)
- The recorded rule that every serde enum is kebab-case: [kb/decisions/cli-config-enum-serialization-kebab-case.md](../decisions/cli-config-enum-serialization-kebab-case.md)

## Impact

- A user cannot copy a command-line value into a config file
- A client that builds the command line from a config, or a config from command-line flags, has to translate the values

## Fix direction

Add `#[serde(rename_all = "kebab-case")]` to the enums the macro generates, regenerate the schemas and the TypeScript client, and check the example configs in `data/` and the pipeline configs for PascalCase values.
