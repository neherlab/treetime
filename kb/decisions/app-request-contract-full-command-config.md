# App clients send the full command config

The HTTP server and the N-API addon accept, for each of the six app commands (`timetree`, `optimize`, `prune`, `ancestral`, `clock`, `mugration`), the same config object that `treetime <command> --config` reads. Every setting the CLI offers is available to the web and desktop apps, and all three clients run a command through one runner in `packages/app-commands`.

## Behavior

- **Request**: the request body is the command's config as JSON, with the keys and values of `packages/schemas/input-config-<command>.schema.json`. The OpenAPI document `packages/app-contracts/openapi.json` embeds these schemas as `<Command>Config` components, and the TypeScript client and zod schemas are generated from it
- **Validation**: the config goes through the CLI config loader: defaults are filled in, the schema check rejects unknown keys (nested ones included), wrong types and invalid enum values, and the conversion to run arguments rejects missing required inputs. Error messages are the ones the CLI prints for the same config file
- **Config check**: `check-config` takes a command name and a config as YAML or JSON text and returns either the config with every default filled in or the error message, its cause chain, and each problem with its location in the text
- **Run events**: every command runs as a run, a record on disk with its own folder. A run emits `started`, then `progress`, `log` and iteration-metric events, then exactly one `terminal` event, and every event is appended to the run's `events.jsonl`. The terminal event is `ok` with the output files that exist after the run and their kinds, `error` with the message and its cause chain, `cancelled`, or `interrupted` for a run that was still running when the process stopped. A rejected config, a failed command and a panic all end in an `error` terminal event. Clients subscribe from an event offset, so a reloaded page or a reopened app resumes a running run's stream
- **Cancellation**: each running run owns its cancellation token, keyed by run id. The server cancels a run on `POST /api/runs/{id}/cancel`; the desktop addon cancels through `RunService.cancel(id)`. Closing an event stream does not cancel the run. Cancelling one run leaves other runs running
- **Paths on the server**: input paths must resolve, after following symbolic links, inside the data directory or the `inputs/` folder of a run (relative paths are resolved against the data directory). Browser files are uploaded into the run's own `inputs/` folder before the run starts. Output paths in the request are discarded; the server writes the outputs of a run to `<runs-dir>/<run-id>/out/`. Translation templates are expanded for each CDS and every expanded path is checked. The addon reads any local path, because the desktop user owns the machine, and writes outputs into the run folder under the app's user-data directory

## Alternatives considered

- **A subset of settings per operation**, with request structs separate from the CLI config. Rejected: the app offers every CLI setting, and separate request structs repeat the fields and defaults of the CLI config and drift from them. The concern that `serde(flatten)` disables `deny_unknown_fields` does not apply, because the schema check runs before deserialization and rejects unknown keys at every level

## Implementation

- `packages/app-commands/src/command.rs`: command names, config preparation, `check-config`, and the runner that returns the written output files
- `packages/app-commands/src/job.rs`: job events, the terminal event and cancellation tokens
- `packages/app-commands/src/runs/`: run records, the run folder layout, event log, lifecycle and output files
- `packages/app-commands/src/config/properties.rs`: path roles (`x-path`) and CLI flags (`x-cli-flag`) of the config settings, read from the schema
- `packages/app-server/src/routes.rs`, `packages/app-server/src/events.rs`, `packages/app-server/src/confine.rs`: HTTP routes, the event stream, and input path confinement
- `packages/app-napi/src/exports.rs`, `packages/app-napi/src/runs.rs`: addon exports and run execution
- `packages/app-contracts/src/bridge.ts`: the TypeScript bridge that turns the terminal event into a result, a `CommandError`, or a `CancelledError`, and exposes runs
