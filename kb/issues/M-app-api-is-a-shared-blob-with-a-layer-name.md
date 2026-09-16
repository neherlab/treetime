# The shared tier is one blob crate named for its layer (`app-api`)

The layer between the core crate `treetime` and the adapters (`app-cli`, `app-server`, `app-napi`) is a single crate, `app-api`, that holds every shared piece at once. It also names itself for its layer position rather than for anything it contains. This violates the shared-tier principles: `kb/PRINCIPLES.md` P3.4 (shared is a tier of cohesive, concern-named units) and P3.7 (no shared blob).

## Symptom

`app-api` mixes unrelated shared concerns in one crate:

- per-command run orchestration and validated argument types (`packages/app-api/src/commands/<cmd>/run.rs`, `.../args.rs`)
- output-path policy and selection (`packages/app-api/src/commands/shared/output.rs`, `.../resolve_outputs.rs`)
- format encoders and writers: tree formats, augur node data, coalescent, GTR, FASTA, CSV, Newick/Nexus comments (`packages/app-api/src/commands/shared/tree_output.rs`, `.../mutation_comment.rs`, and per-command `tree_output.rs` / `augur_node_data.rs`)

The name `app-api` describes none of this; it is a layer/role label.

## Impact

- Any change to one shared concern recompiles the whole tier and forces every adapter to rebuild.
- Ownership is opaque: a reader cannot tell from the crate name or its top level which concern owns a file.
- A single catch-all crate is a magnet for further grab-bag growth.
- Consumers depend on the whole tier rather than only the units they use.

## Fix approach

Split `app-api` into cohesive, concern-named units and retire the `app-api` name. Candidate boundaries:

- tree/format output encoding (per format or one output-encoding unit)
- command configuration: raw + validated argument types and the args-to-core-`Params` mapping
- output-path policy and selection

Preserve the one-way dependency: core (`treetime`) <- shared units <- adapters, and keep adapters independent of each other. Behavior must stay bit-identical.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md): its open design question (a shared application crate vs each adapter direct) is now settled in favor of a shared tier; this issue tracks making that tier well-formed. Its evidence referencing `packages/treetime/src/commands/` is stale once the command layer moves into the shared tier.
- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
- [M-command-output-ownership-is-scattered.md](M-command-output-ownership-is-scattered.md)
- [M-mugration-analysis-interface-exposes-policy-wiring.md](M-mugration-analysis-interface-exposes-policy-wiring.md)
- [M-output-module-mixes-topology-ordering.md](M-output-module-mixes-topology-ordering.md)
