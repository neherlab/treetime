# deser derives slow down rebuilds after an edit

The deser derives generate more code per type than the serde derives did: besides serialization and deserialization, every derived type gets a `describe()` implementation. A rebuild after an edit that reaches `app-commands`, the crate with most derived types and their uses, takes 4 to 5 s longer than with serde.

## Evidence

Build of the `treetime` binary without the compiler cache, the same workspace with serde ("before") and with deser:

| Build                                                                         | Before         | deser          |
| ----------------------------------------------------------------------------- | -------------- | -------------- |
| Clean, release profile (2 runs)                                               | 100 s, 93 s    | 102 s, 96 s    |
| Clean, dev profile                                                            | 150 s          | 166 s          |
| New derived type in `packages/treetime-io/src/auspice_types.rs`, release      | 14.9 to 15.7 s | 20.0 to 20.7 s |
| New derived type in `packages/app-commands/src/results/homoplasy.rs`, release | 13.8 to 15.7 s | 17.3 to 19.5 s |
| New derived type in `packages/treetime-io/src/auspice_types.rs`, dev          | 7.3 to 7.9 s   | 8.6 to 11.5 s  |

The cargo timing report of the `treetime-io` edit in the release profile: `app-commands` takes 16.3 s instead of 11.8 s (code generation 12.0 s instead of 7.2 s), and `treetime-io` 2.1 s instead of 0.6 s. The `treetime-io` library is 7.7 MB instead of 4.0 MB.

> [!IMPORTANT]
> **Investigation required.** Find which generated items cost the code generation time in `app-commands` (the `describe()` implementations, the field visitors, or the instantiation of generic deser code), for example by counting LLVM IR lines per function, and whether deser can skip the parts that TreeTime does not use.
