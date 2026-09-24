# App clients send the full command config

The HTTP server and the N-API addon accept, for each of the six app commands (`timetree`, `optimize`, `prune`, `ancestral`, `clock`, `mugration`), the same config object that `treetime <command> --config` reads. Every setting the CLI offers is available to the web and desktop apps, and all three clients run a command through one runner in `packages/app-commands`.

## Behavior

- **Request**: the request body is the command's config as JSON, with the keys and values of `packages/schemas/input-config-<command>.schema.json`. The OpenAPI document `packages/app-contracts/openapi.json` embeds these schemas as `<Command>Config` components, and the TypeScript client and zod schemas are generated from it
- **Validation**: the config goes through the CLI config loader: defaults are filled in, the schema check rejects unknown keys (nested ones included), wrong types and invalid enum values, and the conversion to run arguments rejects missing required inputs. Error messages are the ones the CLI prints for the same config file
- **Config check**: `check-config` takes a command name and a config as YAML or JSON text and returns either the config with every default filled in or the error message, its cause chain, and each problem with its location in the text
- **Job events**: a run emits `started` (with the job id), then `progress` and `log` events, then exactly one `terminal` event. The terminal event is `ok` with the paths of the output-plan files that exist after the run (the `clock` chart images are not part of the plan), `error` with the message and its cause chain, or `cancelled`. A rejected config, a failed command and a panic all end in an `error` terminal event
- **Cancellation**: each job owns its cancellation token in a registry keyed by job id. The server cancels a job on `POST /api/jobs/{job_id}/cancel` or when the client closes the event stream; the desktop addon cancels through `CommandRunner.cancel(jobId)`. Cancelling one job leaves other jobs running
- **Paths on the server**: input paths must resolve, after following symbolic links, inside the data directory (relative paths are resolved against it). Output paths in the request are discarded; the server writes the outputs of a job to `<out-dir>/<job-id>/`. Translation templates are expanded for each CDS and every expanded path is checked. The addon reads and writes any local path, because the desktop user owns the machine

## Alternatives considered

- **A subset of settings per operation**, with request structs separate from the CLI config. Rejected: the app offers every CLI setting, and separate request structs repeat the fields and defaults of the CLI config and drift from them. The concern that `serde(flatten)` disables `deny_unknown_fields` does not apply, because the schema check runs before deserialization and rejects unknown keys at every level

## Implementation

- `packages/app-commands/src/command.rs`: command names, config preparation, `check-config`, and the runner that returns the written output files
- `packages/app-commands/src/job.rs`: job ids, job events, the terminal event, the job registry and cancellation tokens
- `packages/app-commands/src/config/properties.rs`: path roles (`x-path`) and CLI flags (`x-cli-flag`) of the config settings, read from the schema
- `packages/app-server/src/routes.rs`, `packages/app-server/src/sse.rs`, `packages/app-server/src/confine.rs`: HTTP routes, the event stream, and input path confinement
- `packages/app-napi/src/exports.rs`, `packages/app-napi/src/jobs.rs`: addon exports and job execution
- `packages/app-contracts/src/bridge.ts`: the TypeScript bridge that turns the terminal event into a result, a `CommandError`, or a `CancelledError`
