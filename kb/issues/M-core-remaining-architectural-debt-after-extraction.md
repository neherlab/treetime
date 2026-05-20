# Remaining architectural debt after domain module extraction

After extracting `ancestral/`, `optimize/`, `clock/` from `commands/`, and breaking the `partition/ <-> optimize/` cycle, several structural issues remain. Identified by 6 independent peer reviews (3 Claude, 3 Codex).

## ~~Dependency cycles~~ Resolved

All three dependency cycles have been broken:

### ~~partition/ <-> gtr/~~ Resolved

GTR mutation counting (`infer_gtr_fitch`, `get_mutation_counts_fitch`) moved from `gtr/infer_gtr/` to `ancestral/gtr_inference`. The pure math solver stays in `gtr/infer_gtr/common`. `gtr/` no longer imports from `partition/`.

### ~~partition/ <-> ancestral/~~ Resolved

`resolve_indels_backward/forward` and `Deletion` moved from `ancestral/fitch_indel` to `seq/indel` (pure sequence operations). `PartitionFitch::compress()` convenience method removed; `create_fitch_partition()` in `ancestral/fitch` replaces it. `partition/` no longer imports from `ancestral/`.

### ~~clock/ -> commands/~~ Resolved

`From<BranchSplitArgs>` moved to `clock/run.rs` as `branch_split_to_params` in prior work.

## ~~representation/ renamed to partition/, restructured~~ Resolved

Module renamed from `representation/` to `partition/`, inner `partition/` subdirectory flattened. Topology cleanup (`collapse`, `merge_shared_mutations`, `polytomy_nodes`) moved to `optimize/topology/`. Reroot types and functions moved to `treetime_graph::reroot`. `partition/algo/` retains only `infer_dense`.

## Domain logic trapped in commands/

### ~~commands/timetree/~~ Resolved

Extracted to top-level `timetree/` domain module.

### mugration/ (resolved)

Domain algorithms extracted: GTR refinement moved to `gtr/refinement.rs`, discrete marginal reconstruction converged into shared `partition/marginal_core.rs` + `partition/marginal_discrete.rs`. The old `discrete_marginal.rs`, `gtr_refinement.rs`, `partition/discrete_states.rs`, and `partition/discrete_states.rs` are deleted.

## ~~Visibility and coupling issues~~ Resolved

### ~~28 pub(crate) items in optimize/~~ Resolved

Internal helpers narrowed to `pub(super)` or private. 7 submodules narrowed from `pub` to `pub(super)`. API functions already `pub`. No `pub(crate)` remains in `optimize/`.

### ~~clap derives in domain types~~ Closed

Decision: keep `clap` as a direct dependency. Domain types derive `clap::ValueEnum`/`clap::Args` directly.

## Dead code and cosmetic

- ~~`PartitionMarginal` dead marker trait~~ Resolved
- `hacks/` module contains one 14-line function (`fix_branch_length`). Module name normalizes technical debt. Relocate to `seq/` or inline at 3 call sites.
- `discrete_states.rs` in `partition/` used exclusively by mugration. Move to `mugration/` or future `mugration/` domain module.
