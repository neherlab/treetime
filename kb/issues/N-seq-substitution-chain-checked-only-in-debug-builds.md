# Substitution chains are checked only in debug builds

`compose_substitutions()` checks with `debug_assert_eq!(ps.qry(), cs.reff(), ...)` that the query state of the parent substitution equals the reference state of the child substitution [packages/treetime/src/seq/mutation.rs#L308-L315](../../packages/treetime/src/seq/mutation.rs#L308-L315). Release builds strip the check, so a broken chain silently produces wrong mutation annotations.

## Fix

Replace the `debug_assert_eq!` with a check that runs in all builds and returns an error. `compose_substitutions()` [packages/treetime/src/seq/mutation.rs#L280](../../packages/treetime/src/seq/mutation.rs#L280) already returns `Result`, so the error propagates without an API change.

## Validation

A test composes two substitutions whose states do not chain and expects the error in a release build.
