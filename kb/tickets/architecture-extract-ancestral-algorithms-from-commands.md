# Extract ancestral reconstruction algorithms from commands/ to domain module

## Description

The ancestral reconstruction algorithms in `commands/ancestral/` are domain logic consumed by 7 modules across the codebase. They should be extracted to a domain module outside `commands/`.

## What to move

- `packages/treetime/src/commands/ancestral/fitch.rs` (640 lines) -- core Fitch parsimony algorithm
- `packages/treetime/src/commands/ancestral/marginal.rs` (157 lines) -- marginal reconstruction

Total: 797 lines of domain algorithm code.

## Consumers (7 modules)

1. `commands/ancestral/` -- the command itself (trivially rewired)
2. `commands/optimize/` -- 2 import sites: `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()`
3. `commands/timetree/` -- 6 import sites: `compress_sequences()` (`fitch.rs#L526`), `get_common_length()` (`fitch.rs#L616`), `initialize_marginal()` (`marginal.rs#L21`), `update_marginal()` (`marginal.rs#L41`)
4. `commands/prune/` -- consumer of ancestral reconstruction
5. `representation/` -- reverse dependency (core importing from commands):
   - `representation/partition/fitch_config.rs#L2` -> `get_common_length()` from `commands/ancestral/fitch`
   - `representation/partition/likelihood.rs#L2` -> `get_common_length()` from `commands/ancestral/fitch`
6. `gtr/` -- imports ancestral reconstruction functions
7. `test_utils.rs#L2-L3` -> `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()` from `commands/ancestral/`

## Cross-command import sites for ancestral

- `timetree` -> `ancestral` (6 sites): `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()`
- `optimize` -> `ancestral` (2 sites): `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()`
- `homoplasy` -> `ancestral`: `TreetimeAncestralArgs`

## Target location

Extract to a new domain module outside `commands/` (e.g. `src/ancestral/` or a new crate). The command module should become a thin wrapper: parse CLI args, call domain functions, write output.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
- [N-topology-cleanup-move-merge-shared-mutations.md](../issues/N-topology-cleanup-move-merge-shared-mutations.md) -- topology operation misplacement (related subset)
