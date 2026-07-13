# Optimize pipeline: non-time feature parity with timetree

The `optimize` command is timetree-without-time: branch-length optimization and marginal ancestral reconstruction without dates, clock, or time inference. Several timetree pipeline steps are meaningful without time information but absent from optimize.

## Pipeline comparison

Steps shared by both commands:

| Step                          | optimize                                 | timetree                                             |
| ----------------------------- | ---------------------------------------- | ---------------------------------------------------- |
| Create marginal partition     | `optimize/pipeline.rs`                   | `timetree/pipeline.rs`                               |
| Initialize + update marginal  | `initialize_marginal`, `update_marginal` | same                                                 |
| Normalize GTR rates           | `normalize_partition_rates`              | inside partition init                                |
| ML branch-length optimization | `run_optimize_loop`                      | `optimize_branch_lengths_pre_step` + refinement loop |
| Final marginal pass           | `update_marginal` (post-loop)            | `update_marginal` (post-loop)                        |

Timetree-only steps requiring time (not applicable): date constraints, clock regression, clock filter, timetree inference, coalescent models, relaxed clock, polytomy resolution, confidence intervals.

## Work items

### 1. Generic reroot scoring architecture

Extract the root search algorithm from `clock/find_best_root/` into a new `reroot/` module. Make scoring pluggable via a `RootStats` trait so root-finding works without date information. Prerequisite for optimize reroot support.

- Proposal: [reroot-generic-scoring-architecture.md](reroot-generic-scoring-architecture.md)
- Issues: [kb/issues/M-clock-mindev-wrong-objective.md](../issues/M-clock-mindev-wrong-objective.md) (fix during clock migration)
- Status: implemented. The generic `reroot/` module (`RootStats` trait, `EdgeCostFn<S>`, generic search) and `DivStats` scoring live in `packages/treetime/src/reroot/`.
- Remaining ticket:
  - [reroot-migrate-clock-to-generic-search.md](../tickets/reroot-migrate-clock-to-generic-search.md) -- `ClockStats` wrapping `ClockSet`, update callers, fix MinDev objective, delete `clock/find_best_root/`

### 2. Optimize reroot support

Wire `--reroot`, `--reroot-tips`, and `--keep-root` into optimize. Two-phase BL optimization pattern (optimize -> reroot -> re-optimize). MinDev and tip-based methods only (date-dependent methods rejected at CLI).

- Proposal: [optimize-reroot-support.md](optimize-reroot-support.md)
- Status: implemented. `--reroot`, `--reroot-tips`, and `--keep-root` are wired into optimize with the two-phase branch-length optimization pattern.

### 3. Short branch pruning

Prune internal branches too short to be resolved by the alignment after the optimization loop converges. v0 criterion: `bl < 0.1 * one_mutation` AND `P(zero) > 0.1`.

- Proposal: [optimize-short-branch-pruning.md](optimize-short-branch-pruning.md)
- Tickets:
  - [optimize-add-short-branch-pruning.md](../tickets/optimize-add-short-branch-pruning.md) -- pruning logic, topology collapse, partition reconciliation

### 4. GTR re-estimation in optimization loop

Alternate model inference and BL optimization per ECM iteration. Already proposed as P9, not duplicated here.

- Proposal: P9 in [optimize-convergence-and-robustness.md](optimize-convergence-and-robustness.md)

### 5. Output file control

Specify output tree file paths, select output formats. Orthogonal to the pipeline items above; addresses a separate user request for workflow integration.

- Proposal: [output-format-selection.md](output-format-selection.md) (existing)
- Ticket: implemented on `feat/cli-output-args`

## Execution order

Items 1-2 are sequential (2 depends on 1). Items 3, 4, 5 are independent of each other and of 1-2.

## Related

- [kb/features/optimize.md](../features/optimize.md) -- feature inventory
