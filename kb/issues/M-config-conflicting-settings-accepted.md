# Configs accept settings that the command line rejects as conflicting

The command line rejects conflicting flags through clap (`conflicts_with`, fixed value counts), but a config file, a pipeline step, a server request and a desktop request do not go through clap. `check-config` and a run therefore accept these combinations and the command picks one of the settings without an error.

## Evidence

- Coalescent priors: `coalescent_opt` conflicts with `coalescent` and `coalescent_skyline`, and `coalescent_skyline` with `coalescent` and `coalescent_opt` [packages/app-commands/src/commands/timetree/args.rs#L325](../../packages/app-commands/src/commands/timetree/args.rs#L325), [packages/app-commands/src/commands/timetree/args.rs#L336](../../packages/app-commands/src/commands/timetree/args.rs#L336)
- Root: `keep_root` conflicts with `reroot` and `reroot_tips` in timetree, clock and optimize [packages/app-commands/src/commands/timetree/args.rs#L414](../../packages/app-commands/src/commands/timetree/args.rs#L414), [packages/app-commands/src/commands/clock/args.rs#L210](../../packages/app-commands/src/commands/clock/args.rs#L210), [packages/app-commands/src/commands/optimize/args.rs#L253](../../packages/app-commands/src/commands/optimize/args.rs#L253); optimize `reroot` and `reroot_tips` conflict [packages/app-commands/src/commands/optimize/args.rs#L241](../../packages/app-commands/src/commands/optimize/args.rs#L241)
- Relaxed clock: `--relax` takes exactly two values (`num_args = 2`), while a config accepts `relax: [1.0]` because the schema declares a list of any length
- Reproduction: `treetime timetree --config c.yaml` with `coalescent: 1.0`, `coalescent_opt: true` and `coalescent_skyline: true` (or `keep_root: true` with `reroot: min-dev`, or `relax: [1.0]`) passes the config check and fails only later, when it opens the missing tree file
- The config check validates the JSON schema and the raw-to-run-argument resolution only [packages/app-commands/src/config/load.rs](../../packages/app-commands/src/config/load.rs)

## Impact

- A config that the equivalent command line rejects runs with a setting silently dropped
- The new-analysis form lists the checks `check-config` classifies (`packages/app-commands/src/run_checks.rs`), which has no conflict rule, so it cannot block these runs

## Potential solutions

- O1. Check clap's conflicts and value counts against the settings that differ from their defaults, in the shared config check, so every client reports the CLI's error
- O2. Encode the conflicts as schema constraints (for example `not` over pairs of non-default values), checked by the existing schema check, with a test that they match clap's
- O3. Reject the combinations in the raw-to-run-argument resolution of each command
