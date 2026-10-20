# Timetree does not resolve polytomies by default, and `--keep-polytomies` has no effect

v0 resolves polytomies unless `--keep-polytomies` is given. v1 resolves them only with the opt-in `--resolve-polytomies`, and ignores `--keep-polytomies`. No decision records the changed default.

> [!IMPORTANT]
> **Decision required.** Choose the default polytomy policy of `timetree`. The only v1 method is stochastic, so the default interacts with the seed policy.
>
> - Resolve by default, as v0 does: default output then depends on the run seed, which `SeedArgs::resolve()` draws and logs when `--seed` is absent, whereas v0's greedy default is deterministic
> - Keep polytomies by default: record the divergence from v0 in `kb/decisions/`
>
> Evidence: v0 sets `resolve_polytomies=(not params.keep_polytomies)` [packages/legacy/treetime/treetime/wrappers.py#L485](../../packages/legacy/treetime/treetime/wrappers.py#L485); [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md) approves the method but not the default.

## v0 behavior

- `--keep-polytomies` defaults to False [packages/legacy/treetime/treetime/argument_parser.py#L266-L271](../../packages/legacy/treetime/treetime/argument_parser.py#L266-L271), so v0 resolves polytomies by default
- The default method is greedy and deterministic: `--stochastic-resolve` defaults to False, and `--greedy-resolve` is the default

## v1 behavior

- `--keep-polytomies` and `--resolve-polytomies` are two independent bool flags [packages/app-commands/src/commands/timetree/args.rs#L207-L213](../../packages/app-commands/src/commands/timetree/args.rs#L207-L213). They have no `conflicts_with`, so both can be set at once
- Only `resolve_polytomies` reaches the core and the seed resolution [packages/app-commands/src/commands/timetree/run.rs#L68](../../packages/app-commands/src/commands/timetree/run.rs#L68) [packages/app-commands/src/commands/timetree/run.rs#L150](../../packages/app-commands/src/commands/timetree/run.rs#L150); the round skips resolution unless it is set [packages/treetime/src/timetree/round.rs#L127-L129](../../packages/treetime/src/timetree/round.rs#L127-L129)
- `keep_polytomies` is never read. The app shows it as the switch "Keep polytomies unresolved" [packages/app-ui/src/analysis/MainSettings.tsx#L213](../../packages/app-ui/src/analysis/MainSettings.tsx#L213) and passes it to the back end [packages/app-backend/src/command.rs#L111](../../packages/app-backend/src/command.rs#L111), so the switch has no effect
- The only method is stochastic ([kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md))

## Impact

- A default `treetime timetree` run keeps multifurcations that v0 resolves, which changes node dates and the coalescent estimates
- A CLI flag and an app switch silently do nothing

## Fix

After the decision, represent polytomy handling as one exhaustive policy enum and parse the CLI flags into it. When both flags are given, clap and every other transport (config files, server requests) reject the request with the same explanation, so the pipeline never receives contradictory booleans.

## Validation

- A default run follows the decided policy
- CLI parsing rejects both flags together and accepts each alone; config files and server requests enforce the same rule
- An exhaustive match prevents a new policy variant from being ignored
- The app switch changes the output

## Related issues

- [M-cli-flags-parsed-but-ignored.md](M-cli-flags-parsed-but-ignored.md)
- [M-timetree-coalescent-modes-represent-invalid-states.md](M-timetree-coalescent-modes-represent-invalid-states.md)
- [M-timetree-polytomy-cleanup-drops-dated-single-child-node.md](M-timetree-polytomy-cleanup-drops-dated-single-child-node.md)
- [M-timetree-relaxed-clock-rates-reset-after-polytomy-resolution.md](M-timetree-relaxed-clock-rates-reset-after-polytomy-resolution.md)
- [H-timetree-polytomy-resolution-blows-up-the-time-grid.md](H-timetree-polytomy-resolution-blows-up-the-time-grid.md)
