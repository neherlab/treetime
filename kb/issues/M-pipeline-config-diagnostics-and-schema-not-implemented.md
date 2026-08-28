# Pipeline and per-command config lack rich diagnostics and an editor schema

## Symptom and reproduction

The `treetime pipeline` command loads and runs a JSON/YAML pipeline config, and every command accepts a `--config` file. Both report configuration errors as plain single-line messages. There is no source snippet, no caret pointing at the offending key, and no batching: a config with several independent mistakes surfaces only the first one, so fixing it takes several iterations.

There is also no generated JSON Schema for either the pipeline config or the per-command config, so an editor offers no completion, hover documentation, or inline validation while writing a config file.

Reproduce by running `treetime pipeline --config <file>` with, for example, a bad enum value and an unknown field in the same step: only one error prints, as text, with no location.

## Impact and scope

- Config authoring is harder than intended: errors are accurate but give no location and no did-you-mean beyond command tags and chained references.
- No editor assistance for pipeline or per-command config files.
- The per-command `--config` overlay path is unchanged and does not share any diagnostics layer, so any future diagnostics improvement must be wired into both entry points.

Functionality is unaffected: pipelines resolve, chain, validate, and run correctly, and configuration errors are reported (just without source carets, batching, or schema).

## Root cause

The configuration diagnostics and schema layers are not implemented. The loader deserializes the config into a permissive value and resolves it, returning the first error as an `eyre` report. No JSON Schema is generated from the command argument types, and no batched, source-annotated rendering is performed.

## Fix approach

Three connected pieces, in order:

1. Derive `JsonSchema` across the command argument structs, the shared argument groups (alignment, alphabet, model, output, gap-fill, metadata, reroot, topology-order, config), and the scientific parameter enums and structs they reference (for example the substitution-model, alphabet, gap-fill, ancestral, optimize, and clock parameter types). Generate `pipeline.schema.json` and a strict `<command>.schema.json` per command from the `schema` subcommand in `app-cli`, applying a schema transform that lets a whole-value template string (`"{{ ... }}"`) satisfy an otherwise scalar field in the pipeline schema only.
2. Build one shared diagnostics module that reads the raw config text, validates the parsed value against the generated schema to collect all structural errors in one pass, builds a byte-span index from the source, validates interpolation references, adds edit-distance suggestions, and renders every collected error as a batched, caret-annotated report to stderr before returning a short summary error.
3. Refactor the per-command `--config` overlay onto the shared diagnostics module so a bad per-command config renders the same caret errors, keeping the precedence (command-line flag over file over default) unchanged.

## Workaround

Use `treetime pipeline --check` to resolve and print the plan (steps, resolved inputs, outputs) without running anything; it surfaces most configuration mistakes with actionable messages. Validate the config against the produced errors iteratively.
