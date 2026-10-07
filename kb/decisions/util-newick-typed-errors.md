# Typed errors in util-newick

Decided on 2026-10-07.

TreeTime reports errors with `eyre`, and the `make_error!` macros of `treetime-utils` build them. `util-newick` is an exception: its writer, its graph model, and its number formatting return the typed error `NewickWriteError`, and parsing a dialect name returns `ParseDialectError`, both derived with `thiserror`. Its reader already returned the typed `NewickError`.

## Reason

`util-newick` has a second user. [TreeKnit](https://github.com/neherlab/treeknit-rs) reads and writes its Newick trees and its ARG output (extended Newick with BEAST annotations) through a copy of the crate. TreeKnit allows `eyre` only in its command line, because its libraries, which also run in the browser as WebAssembly, return typed errors. With typed errors, the two copies stay the same, and a change to the crate can be copied without edits.

## Scope

- `util-newick` only. The other crates of TreeTime keep `eyre`
- `NewickWriteError` holds one message. A step that adds context prefixes it, as `wrap_err` does, so the text reads the same as the `eyre` chain did: "When writing Newick: When writing the branch above node 1 ('A'): Newick cannot represent the number inf"
- TreeTime code converts the errors into `eyre::Report` with `?`, as for any error type
- The crate depends on no other crate of TreeTime, so it can be copied as it is
