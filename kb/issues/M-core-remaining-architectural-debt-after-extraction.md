# Remaining architectural debt after domain module extraction

After extracting `ancestral/`, `optimize/`, `clock/` from `commands/`, and breaking the `representation/ <-> optimize/` cycle, several structural issues remain. Identified by 6 independent peer reviews (3 Claude, 3 Codex).

## Dependency cycles

Three bidirectional dependency cycles remain between top-level modules:

### representation/ <-> gtr/ (10+ import lines)

Every partition file imports `GTR` from `gtr/`. Two files in `gtr/infer_gtr/` import back from `representation/`:

- `gtr/infer_gtr/dense.rs` -> `representation::partition::marginal_dense::PartitionMarginalDense`
- `gtr/infer_gtr/fitch.rs` -> `representation::partition::traits::PartitionCompressed`

Breakable: parameterize GTR inference over iterators instead of concrete partition types. `infer_gtr/dense.rs` reads `(msg_to_parent, msg_to_child)` pairs - accepting `impl Iterator<Item = (&Array2<f64>, &Array2<f64>)>` breaks the reverse direction.

### representation/ <-> ancestral/ (2 upward imports)

- `partition/fitch.rs:2:` -> `ancestral::fitch::compress_sequences`
- `partition/marginal_dense.rs:4:` -> `ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward}`

Breakable: `resolve_indels_*` are pure functions on gap ranges (`Vec<(usize, usize)>`, `BTreeMap<(usize, usize), Deletion>`) - move to `seq/` or `representation/`. `compress_sequences` is called by `PartitionFitch::compress()` convenience constructor - remove the convenience method or invert the call direction.

### clock/ -> commands/ (1 import)

- `clock/find_best_root/params.rs:1:` -> `commands::clock::args::BranchSplitArgs`

The `From<&BranchSplitArgs> for BranchPointOptimizationParams` impl lives in the domain module. Move it to `commands/clock/`.

## Domain logic trapped in commands/

### ~~commands/timetree/~~ **Resolved**

Extracted to top-level `timetree/` domain module: `inference/` (backward/forward passes, branch length likelihood, runner), `optimization/` (polytomy, relaxed clock, clock filter, reroot), `utils.rs`. `commands/timetree/` retains CLI-coupled modules (args, run, output, convergence, initialization, refinement).

### commands/mugration/ (1227 lines)

Contains domain algorithms:

- `gtr_refinement.rs` (381 lines) - iterative GTR inference from discrete trait data
- `discrete_marginal.rs` (170 lines) - discrete-trait marginal reconstruction

These import only from `gtr/` and `partition/discrete` - no args coupling.

## representation/ is not a coherent module

5100 lines mixing three concerns under a vague name:

- `payload/` (1385 lines) - node/edge data types. Pure data.
- `partition/` (2956 lines) - partition types AND marginal inference algorithms (1811 lines of science code) AND optimization coefficient types (268 lines)
- `algo/` (603 lines) - topology cleanup operations

The marginal inference algorithms (`marginal_dense.rs` 711 lines, `marginal_passes.rs` 432 lines, `marginal_helpers.rs` 316 lines) are the densest scientific code in the project. They live in `representation/` as trait implementations on partition types while `ancestral/marginal.rs` (160 lines) is a thin tree-walking orchestrator. The architecture is inverted: the "domain module" wraps, the "data module" computes.

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
- `clock/reroot.rs` re-exports `EdgeMergeInfo`, `EdgeSplitInfo`, `RerootChanges`, `RerootResult` from `partition/algo/topology_cleanup/reroot`. Creates dual import paths. No production caller uses the re-export. Remove.
- `hacks/` module contains one 14-line function (`fix_branch_length`). Module name normalizes technical debt. Relocate to `seq/` or inline at 3 call sites.
- `discrete_states.rs` in `representation/` used exclusively by mugration. Move to `commands/mugration/` or future `mugration/` domain module.
