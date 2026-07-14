# Library crate is not consumable by non-CLI clients

The `treetime` library crate cannot be used by Python bindings, desktop apps, or web frontends. The library depends on `clap` (CLI framework), contains CLI argument structs in domain modules, and has domain logic trapped inside `commands/` (a CLI-layer concept). Any non-CLI consumer inherits the full CLI dependency graph and must construct clap-derived argument structs to call domain functions.

This is a structural blocker for multi-client architecture: `treetime-python` (PyO3 bindings), `treetime-desktop` (Electron/Tauri wrapper), or any future consumer.

## Evidence

`commands/` in the library crate contains approximately 15,000 lines that mix domain algorithms with CLI orchestration, and domain types derive `clap::Args` or `clap::ValueEnum`. The boundary therefore makes CLI-layer types and dependencies prerequisites for non-CLI consumers.

## Scope

### ~~clap contamination in library domain types~~ Closed

Decision: keep `clap` as a direct dependency of the `treetime` library crate. Domain types derive `clap::ValueEnum`/`clap::Args` directly. Non-CLI clients accept the transitive dependency.

### Domain logic in commands/

`commands/timetree/` (~10k lines) and `mugration/` (~2.2k lines) contain domain algorithms that should be top-level modules. Tracked in detail by:

- `H-core-command-module-shared-ops-entanglement.md` (remaining items)
- `M-core-remaining-architectural-debt-after-extraction.md` (timetree, mugration, representation)

### Missing client crates

- `treetime-python`: PyO3 bindings + Python CLI wrapping Rust library (not started)
- ~~`treetime-desktop`: thin wrapper with minimal desktop frontend~~ (scaffolded: napi-rs async wrappers for all commands, Electron frontend placeholder)

## Impact

- Encourages contributors and automation clients to treat `treetime-cli` as the primary crate, moving domain logic out of the library

## Related issues

- [Command modules contain shared operations that belong in domain layers](H-core-command-module-shared-ops-entanglement.md)
- [Remaining architectural debt after domain module extraction](M-core-remaining-architectural-debt-after-extraction.md)

## Related tickets

- [kb/tickets/architecture-add-treetime-python-crate.md](../tickets/architecture-add-treetime-python-crate.md)
- [kb/tickets/architecture-migrate-command-tests-with-dissolution.md](../tickets/architecture-migrate-command-tests-with-dissolution.md)
- [kb/tickets/architecture-move-commands-to-cli-crate.md](../tickets/architecture-move-commands-to-cli-crate.md)
- [kb/tickets/architecture-restructure-packages-for-multi-client.md](../tickets/architecture-restructure-packages-for-multi-client.md)
