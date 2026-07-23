# Graph capability contracts silently discard state

`TimeLength` combines reading and writing an optional time length [`packages/treetime-graph/src/edge.rs#L38`](../../packages/treetime-graph/src/edge.rs#L38). `EdgeClock` implements that capability by always returning `None` and ignoring every setter call [`packages/treetime/src/clock/clock_graph.rs#L160`](../../packages/treetime/src/clock/clock_graph.rs#L160).

`Divergence` similarly accepts `Option<f64>`, while its own contract permits non-optional implementations to ignore `set_div(None)` [`packages/treetime-graph/src/node.rs#L27`](../../packages/treetime-graph/src/node.rs#L27). A generic caller can therefore issue a supported command and immediately observe that the subtype refused it.

## Contract violation

A caller constrained by `TimeLength` is entitled to expect that `set_time_length(value)` changes the subsequently observed value. `EdgeClock` accepts the command but discards it. `NodeTimetree` similarly has required divergence storage while implementing an optional-state setter whose `None` command cannot be honored.

These no-op overrides hide unsupported state transitions instead of expressing them in the type system. Algorithms can appear generic while silently producing incomplete payloads for one implementation.

## Required trait split

- Read-only access must be independent from writable access.
- Optional storage may accept `Option<T>` setters.
- Required storage must accept `T` and cannot advertise removal.
- Types without a value must not implement the corresponding storage capability.

Callers that genuinely support several storage contracts must match those contracts explicitly rather than relying on a common no-op setter.

## Validation

- Trait-law tests verify every writable implementation observes the value just set.
- Compile-time bounds exclude `EdgeClock` from operations requiring writable time length.
- Required divergence cannot be cleared through any public API.
- Clock, timetree, reroot, and serialization tests prove unchanged domain behavior.
