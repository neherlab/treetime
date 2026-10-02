# Command orchestration mixes application policy, domain workflow, and I/O

The command boundary is not a thin adapter. Command runners combine argument translation, input loading, pipeline invocation, output policy, serialization, and filesystem writes, while one timetree pipeline function sequences the complete scientific workflow.

## Evidence

- `fn run_ancestral_reconstruction()` reads FASTA, maps CLI arguments, runs inference, builds output projections, and writes several formats [packages/app-commands/src/commands/ancestral/run.rs#L42](../../packages/app-commands/src/commands/ancestral/run.rs#L42).
- The same application-level shape appears in clock, mugration, optimize, prune, and timetree runners under [`packages/app-commands/src/commands`](../../packages/app-commands/src/commands).
- `fn timetree::pipeline::run()` [packages/treetime/src/timetree/pipeline.rs#L55-L141](../../packages/treetime/src/timetree/pipeline.rs#L55-L141) sequences the complete scientific workflow through step functions: date loading, clock estimation, the pre-loop steps (ML branch-length optimization, rerooting, clock filter), coalescent initialization, refinement, confidence intervals, and result assembly.
Domain modules do not import `commands/`, but the remaining application orchestration has no explicit owner.

## Open design question

The application operation boundary must serve CLI, HTTP, N-API, desktop, and future Python clients. Moving all orchestration into the CLI crate would leave the other clients without an owner. The unresolved choice is whether a shared application crate owns validated operations and in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types.

No ticket should move `commands/` until this boundary is decided. Scientific stage ordering and fallback behavior must remain unchanged unless separately approved.

## Related issues

- [M-core-partition-init-orchestration-duplication.md](M-core-partition-init-orchestration-duplication.md)
- [M-command-output-ownership-is-scattered.md](M-command-output-ownership-is-scattered.md)
- [M-mugration-analysis-interface-exposes-policy-wiring.md](M-mugration-analysis-interface-exposes-policy-wiring.md)
- [M-output-module-mixes-topology-ordering.md](M-output-module-mixes-topology-ordering.md)
