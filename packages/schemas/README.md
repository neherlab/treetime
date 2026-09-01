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

For a YAML config, add a `yaml-language-server` modeline as the first line, pointing at the matching schema. The editor reads it for completion and validation.

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-pipeline.schema.json
```

The example configs under `data/` carry such a modeline pointing at the `rust` branch. A JSON config, which cannot carry a comment, uses a top-level `$schema` key with the same URL instead; the loader accepts and ignores it. The [Red Hat YAML extension](https://marketplace.visualstudio.com/items?itemName=redhat.vscode-yaml) resolves a relative path from the config file's location if you prefer a local checkout to the branch URL.

## Regenerating

These files are generated from the command argument types. Regenerate them after changing any command's arguments or the pipeline shape:

```bash
treetime schema --for all -o packages/schemas
```

A drift-guard test in `app-cli` fails if a committed schema here no longer matches the code.
