# Remaining architectural debt after domain module extraction

After extracting `ancestral/`, `optimize/`, `clock/` from `commands/`, and breaking the `representation/ <-> optimize/` cycle, several structural issues remain. Identified by 6 independent peer reviews (3 Claude, 3 Codex).

## ~~Dependency cycles~~ Resolved

All three dependency cycles have been broken:

### ~~partition/ <-> gtr/~~ Resolved

GTR mutation counting (`infer_gtr_fitch`, `get_mutation_counts_fitch`) moved from `gtr/infer_gtr/` to `ancestral/gtr_inference`. The pure math solver stays in `gtr/infer_gtr/common`. `gtr/` no longer imports from `partition/`.

### ~~partition/ <-> ancestral/~~ Resolved

`resolve_indels_backward/forward` and `Deletion` moved from `ancestral/fitch_indel` to `seq/indel` (pure sequence operations). `PartitionFitch::compress()` convenience method removed; `create_fitch_partition()` in `ancestral/fitch` replaces it. `partition/` no longer imports from `ancestral/`.

### ~~clock/ -> commands/~~ Resolved

`From<BranchSplitArgs>` moved to `commands/clock/run.rs` as `branch_split_to_params` in prior work.

## ~~representation/ renamed to partition/, restructured~~ Resolved

Module renamed from `representation/` to `partition/`, inner `partition/` subdirectory flattened. Topology cleanup (`collapse`, `merge_shared_mutations`, `polytomy_nodes`) moved to `optimize/topology/`. Reroot types and functions moved to `treetime_graph::reroot`. `partition/algo/` retains only `infer_dense`.

## Domain logic trapped in commands/

### ~~commands/timetree/~~ Resolved

Extracted to top-level `timetree/` domain module.

### commands/mugration/ (1227 lines)

Contains domain algorithms:

- `gtr_refinement.rs` (381 lines) - iterative GTR inference from discrete trait data
- `discrete_marginal.rs` (170 lines) - discrete-trait marginal reconstruction

These import only from `gtr/` and `partition/discrete` - no args coupling.

## Visibility and coupling issues

### 28 pub(crate) items in optimize/

Most are internal helpers (`newton_tolerance_*`, `chain_rule_*`, `brent_bracket`, etc.) that no code outside `optimize/` calls. True API functions (`run_optimize_mixed`, `evaluate_with_indels_log_lh_only`) used by timetree need `pub` for future crate split.

### clap derives in domain types

6 domain-layer files have `use clap::ValueEnum` or `use clap::Args`:

- `clock/clock_regression.rs` - `ClockParams` has `#[derive(Args)]`
- `clock/find_best_root/params.rs` - 5 types with `ValueEnum`/`Args`
- `optimize/args.rs` - `BranchOptMethod`, `InitialGuessMode`
- `alphabet/alphabet.rs` - `AlphabetName`
- `gtr/get_gtr.rs` - `GtrModelName`
- `seq/gap_fill.rs` - `GapFillMode`

Domain types should not depend on CLI parsing. Library consumers pull in `clap` as transitive dependency.

## Dead code and cosmetic

- `PartitionMarginal` in `partition/traits.rs:115:` is a dead marker trait. Empty impl on both partition types, used only as supertrait bound on `PartitionMarginalOps`. Safe to remove.
- `clock/reroot.rs` re-exports `EdgeMergeInfo`, `EdgeSplitInfo`, `RerootChanges`, `RerootResult` from `treetime_graph::reroot`. Creates dual import paths. No production caller uses the re-export. Remove.
- `hacks/` module contains one 14-line function (`fix_branch_length`). Module name normalizes technical debt. Relocate to `seq/` or inline at 3 call sites.
- `discrete_states.rs` in `partition/` used exclusively by mugration. Move to `commands/mugration/` or future `mugration/` domain module.
