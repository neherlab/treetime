# Timetree does not resolve polytomies by default

v0 resolves polytomies unless `--keep-polytomies` is given. v1 resolves them only with the opt-in `--resolve-polytomies`, and ignores `--keep-polytomies`. No decision records the changed default.

## v0 behavior

- `resolve_polytomies=(not params.keep_polytomies)` ([packages/legacy/treetime/treetime/wrappers.py#L485](../../packages/legacy/treetime/treetime/wrappers.py#L485)); `--keep-polytomies` defaults to False ([packages/legacy/treetime/treetime/argument_parser.py#L266-L271](../../packages/legacy/treetime/treetime/argument_parser.py#L266-L271))
- The default method is greedy and deterministic: `--stochastic-resolve` defaults to False, and `--greedy-resolve` is the default

## v1 behavior

- `--resolve-polytomies` is an opt-in bool ([packages/app-cli/src/commands/timetree/args.rs](../../packages/app-cli/src/commands/timetree/args.rs)) read by `refine_topology` in [packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs)
- `keep_polytomies` is copied into `TimetreeParams` and never read
- The only method is stochastic ([kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md), which approves the method but not the default). The seed is random unless `--seed` is given

## Impact

A default `treetime timetree` run keeps multifurcations that v0 resolves, which changes node dates and the coalescent estimates. Resolving by default would make default output depend on a random seed unless a fixed default seed is used, whereas v0's greedy default is deterministic.

## Related issues

- [N-timetree-polytomy-flags-no-conflict.md](N-timetree-polytomy-flags-no-conflict.md)
