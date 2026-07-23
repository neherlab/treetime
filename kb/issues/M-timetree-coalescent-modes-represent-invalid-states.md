# Timetree parameters can represent conflicting coalescent modes

CLI parsing rejects skyline together with constant or optimized coalescent modes [`packages/treetime/src/commands/timetree/args.rs#L121-L140`](../../packages/treetime/src/commands/timetree/args.rs#L121-L140). The transport-neutral `TimetreeParams` represents those modes as `Option<f64>` plus independent booleans, so server or direct callers can construct `coalescent_skyline` together with `coalescent` or `coalescent_opt` [`packages/treetime/src/timetree/pipeline.rs#L65-L67`](../../packages/treetime/src/timetree/pipeline.rs#L65-L67).

`fn coalescent_initialization()` resolves the booleans by precedence instead of receiving a parsed mode [`packages/treetime/src/timetree/pipeline.rs#L444-L454`](../../packages/treetime/src/timetree/pipeline.rs#L444-L454).

## Invalid combinations

The parameter struct can simultaneously request skyline and a constant timescale, skyline and optimization, or several booleans whose precedence determines the result. These states are impossible through one CLI parser path but remain constructible by server conversions, tests, and direct library callers.

Validation confined to clap is therefore insufficient for the transport-neutral operation contract. Precedence in the pipeline silently chooses one interpretation instead of rejecting an invalid request.

## Required type

Represent the modes as one exhaustive value:

- disabled;
- constant with a required timescale;
- optimized with an explicit optional initial value;
- skyline with its skyline-specific configuration.

Parse CLI and transport DTOs into that type before entering the pipeline. Keep the scientific meaning and defaults of each current valid mode unchanged.

## Validation

- Every transport adapter rejects the same invalid source combinations.
- Exhaustive unit cases cover all valid modes without precedence logic.
- Golden-master timetree tests preserve constant, optimized, skyline, and disabled behavior.
- Serialization round trips cannot construct a conflicting internal state.

## Related issues

- [N-timetree-polytomy-flags-no-conflict.md](N-timetree-polytomy-flags-no-conflict.md)
- [H-app-transport-contracts-diverge-across-clients.md](H-app-transport-contracts-diverge-across-clients.md)
