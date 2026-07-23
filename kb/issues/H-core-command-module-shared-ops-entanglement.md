# Command orchestration mixes application policy, domain workflow, and I/O

The command boundary is not a thin adapter. Command runners combine argument translation, input loading, pipeline invocation, output policy, serialization, and filesystem writes, while the timetree pipeline concentrates most of the scientific workflow in one function.

## Evidence

- `fn run_ancestral_reconstruction()` reads FASTA, maps CLI arguments, runs inference, builds output projections, and writes several formats [`packages/treetime/src/commands/ancestral/run.rs#L31`](../../packages/treetime/src/commands/ancestral/run.rs#L31).
- The same application-level shape appears in clock, mugration, optimize, prune, and timetree runners under [`packages/treetime/src/commands`](../../packages/treetime/src/commands).
- `fn timetree::pipeline::run()` spans the complete scientific sequence from date loading and clock estimation through coalescent initialization, refinement, rerooting, confidence intervals, and result assembly [`packages/treetime/src/timetree/pipeline.rs#L102`](../../packages/treetime/src/timetree/pipeline.rs#L102).
- `fn run_refinement_iteration()` combines relaxed-clock application, topology resolution, partition reconciliation, ancestral-state comparison, and clock re-estimation [`packages/treetime/src/timetree/refinement.rs#L27`](../../packages/treetime/src/timetree/refinement.rs#L27).

Domain modules do not import `commands/`, but the remaining application orchestration has no explicit owner and the timetree workflow does not expose its scientifically meaningful state transitions.

## Open design question

The application operation boundary must serve CLI, HTTP, N-API, desktop, and future Python clients. Moving all orchestration into the CLI crate would leave the other clients without an owner. The unresolved choice is whether a shared application crate owns validated operations and in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types.

No ticket should move `commands/` until this boundary is decided. Scientific stage ordering and fallback behavior must remain unchanged unless separately approved.

## Related issues

- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
- [H-app-transport-contracts-diverge-across-clients.md](H-app-transport-contracts-diverge-across-clients.md)
- [M-core-partition-init-orchestration-duplication.md](M-core-partition-init-orchestration-duplication.md)
- [M-command-output-ownership-is-scattered.md](M-command-output-ownership-is-scattered.md)
- [M-mugration-analysis-interface-exposes-policy-wiring.md](M-mugration-analysis-interface-exposes-policy-wiring.md)
- [M-timetree-refinement-iteration-mixes-state-transitions.md](M-timetree-refinement-iteration-mixes-state-transitions.md)
- [M-output-module-mixes-topology-ordering.md](M-output-module-mixes-topology-ordering.md)
