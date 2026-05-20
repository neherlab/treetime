# Library crate is not consumable by non-CLI clients

The `treetime` library crate cannot be used by Python bindings, desktop apps, or web frontends. The library depends on `clap` (CLI framework), contains CLI argument structs in domain modules, and has domain logic trapped inside `commands/` (a CLI-layer concept). Any non-CLI consumer inherits the full CLI dependency graph and must construct clap-derived argument structs to call domain functions.

This is a structural blocker for multi-client architecture: `treetime-python` (PyO3 bindings), `treetime-desktop` (Electron/Tauri wrapper), or any future consumer.

## Evidence

An AI refactoring agent moved the entire timetree implementation from the library into the CLI crate, treating `treetime-cli` as the owner of domain logic. This reflects the current blurred boundary: `commands/` in the library crate contains ~15k lines mixing domain algorithms with CLI orchestration, and domain types derive `clap::Args` / `clap::ValueEnum`.

## Scope

### ~~clap contamination in library domain types~~ Closed

Decision: keep `clap` as a direct dependency of the `treetime` library crate. Domain types derive `clap::ValueEnum`/`clap::Args` directly. Non-CLI clients accept the transitive dependency.

### Domain logic in commands/

`commands/timetree/` (~10k lines) and `mugration/` (~2.2k lines) contain domain algorithms that should be top-level modules. Tracked in detail by:

- `H-core-command-module-shared-ops-entanglement.md` (remaining items)
- `M-core-remaining-architectural-debt-after-extraction.md` (timetree, mugration, representation)

### Missing client crates

No crates exist for non-CLI consumption:

- `treetime-python`: PyO3 bindings + Python CLI wrapping Rust library
- `treetime-desktop`: thin wrapper with minimal desktop frontend

These crates create structural backpressure: if the library must remain consumable by multiple clients, domain logic cannot drift into any single client crate.

## Impact

- Encourages AI agents and contributors to treat `treetime-cli` as the primary crate, moving domain logic out of the library

## Related issues

- [Command modules contain shared operations that belong in domain layers](H-core-command-module-shared-ops-entanglement.md)
- [Remaining architectural debt after domain module extraction](M-core-remaining-architectural-debt-after-extraction.md)
