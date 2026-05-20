# Add treetime-python crate with PyO3 bindings and Python CLI

## Description

Create `packages/treetime-python/` as a workspace crate providing Python bindings via PyO3 and a Python CLI wrapping the Rust library. This crate depends on the `treetime` library crate and exposes domain functions as a Python module.

The crate serves two purposes:

1. Structural backpressure: its existence as a library consumer prevents domain logic from migrating into `treetime-cli`. If a function is needed by both CLI and Python, it must live in the library.
2. Python ecosystem access: replaces the legacy Python implementation (`packages/legacy/treetime/`) with a Rust-backed Python package using modern Python tooling.

## Crate structure

```
packages/treetime-python/
  Cargo.toml          # PyO3 cdylib + workspace inheritance
  pyproject.toml      # maturin build backend, uv/hatch project
  python/
    treetime/
      __init__.py     # re-exports from Rust module
      cli.py          # click/typer CLI wrapping Rust functions
  src/
    lib.rs            # PyO3 module definition
    ancestral.rs      # Python-facing ancestral API
    clock.rs          # Python-facing clock API
    timetree.rs       # Python-facing timetree API
    mugration.rs      # Python-facing mugration API
    types.rs          # Python type conversions (PyGraph, PyGTR, etc.)
```

## Dependencies

- `treetime` (workspace) - core library
- `treetime-io` (workspace) - file format I/O
- `treetime-graph` (workspace) - graph types
- `treetime-utils` (workspace) - error handling
- `pyo3` - Rust-Python bindings
- `numpy` (via `pyo3`) - ndarray <-> numpy conversion

Python side:

- `maturin` - build backend
- `click` or `typer` - Python CLI framework
- `numpy` - array interop

## Workspace integration

- Add to `[workspace.members]` in root `Cargo.toml`
- Add `pyo3` and `numpy` to `[workspace.dependencies]`
- CI: build wheel, run Python tests

## Scope

Initial implementation should expose a minimal API surface covering one command (e.g., `ancestral`) end-to-end. Expand incrementally as the library API stabilizes.

## Prerequisites

- `commands/` module moved to `treetime-cli` (ticket: `architecture-move-commands-to-cli-crate.md`)
- timetree domain logic extracted from `commands/` (ticket: `architecture-extract-timetree-inference-from-commands.md`)

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
