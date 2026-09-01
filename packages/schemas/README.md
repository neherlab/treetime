# JSON schemas

JSON Schema (draft 2020-12) documents for TreeTime's data shapes. Editors use them to give completion, hover documentation, and inline validation while you edit a config file.

## Naming

Files are named `<category>-<name>.schema.json` so the category is legible from the filename and related schemas sort together.

- `input-config-*` -- the argument object of a command run with `--config`, and the `pipeline` config. `<name>` is the command, or `pipeline` for the whole-pipeline flavor.

Further categories (other inputs, and command outputs) will follow the same `<category>-<name>` scheme.

| Schema                                                                     | Applies to                           |
| -------------------------------------------------------------------------- | ------------------------------------ |
| [`input-config-pipeline.schema.json`](input-config-pipeline.schema.json)   | A `treetime pipeline --config` file  |
| [`input-config-ancestral.schema.json`](input-config-ancestral.schema.json) | A `treetime ancestral --config` file |
| [`input-config-clock.schema.json`](input-config-clock.schema.json)         | A `treetime clock --config` file     |
| [`input-config-mugration.schema.json`](input-config-mugration.schema.json) | A `treetime mugration --config` file |
| [`input-config-optimize.schema.json`](input-config-optimize.schema.json)   | A `treetime optimize --config` file  |
| [`input-config-prune.schema.json`](input-config-prune.schema.json)         | A `treetime prune --config` file     |
| [`input-config-timetree.schema.json`](input-config-timetree.schema.json)   | A `treetime timetree --config` file  |

The runtime data-contract schemas (`version-info`, `progress-event`, `error-response`) are also emitted here. Their source of truth for the TypeScript bindings is `packages/app-contracts/src/generated/`; the copies here keep this directory a complete schema set.

## Associating a schema with a config file

Add a top-level `$schema` key whose value is the URL of the matching schema. It works in both YAML and JSON, the editor reads it for completion and validation, and the config loader accepts and ignores it.

```yaml
"$schema": "https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-pipeline.schema.json"
```

The example configs under `data/` carry such a key pointing at the `rust` branch. For a local checkout you may instead point at a relative path (resolved from the config file's location by the [Red Hat YAML extension](https://marketplace.visualstudio.com/items?itemName=redhat.vscode-yaml)).

## Regenerating

These files are generated from the command argument types. Regenerate them after changing any command's arguments or the pipeline shape:

```bash
treetime schema --for all -o packages/schemas
```

A drift-guard test in `app-cli` fails if a committed schema here no longer matches the code.
