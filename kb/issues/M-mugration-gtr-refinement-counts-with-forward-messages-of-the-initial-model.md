# Mugration GTR refinement counts transitions with the forward messages of the initial model

In every iteration of the mugration GTR refinement, v1 counts transitions with the forward messages of the initial uniform model. v0 counts with the forward messages of the first fitted GTR. The fitted model, the rate, and the trait confidences can therefore differ from v0 on default settings.

## v1 behavior

Mugration builds the initial GTR with a uniform exchangeability matrix and runs one marginal update under it [packages/treetime/src/mugration/pipeline.rs#L99-L120](../../packages/treetime/src/mugration/pipeline.rs#L99-L120). `fn refine_gtr_model_and_rate()` takes the forward messages of that update once [packages/treetime/src/gtr/refinement.rs#L73-L77](../../packages/treetime/src/gtr/refinement.rs#L73-L77). In each of its iterations, it counts transitions with the current node states and backward messages and with these initial forward messages, fits a new GTR, and optimizes the rate [packages/treetime/src/gtr/refinement.rs#L87-L99](../../packages/treetime/src/gtr/refinement.rs#L87-L99). `fn optimize_gtr_rate()` runs only backward passes, so it returns new node states and backward messages but no forward messages. The dense count reads both message directions [packages/treetime/src/partition/marginal/shared/data.rs#L33-L37](../../packages/treetime/src/partition/marginal/shared/data.rs#L33-L37).

## v0 behavior

`reconstruct_discrete_traits()` calls `infer_ancestral_sequences(infer_gtr=True, marginal=True)` and then `optimize_gtr_rate()` once, and then repeats `infer_gtr()` and `optimize_gtr_rate()` for each iteration [packages/legacy/treetime/treetime/wrappers.py#L785-L802](../../packages/legacy/treetime/treetime/wrappers.py#L785-L802):

1. `infer_ancestral_sequences()` first fits a GTR from a marginal reconstruction under the initial model, then runs a full marginal reconstruction under that first fitted GTR [packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566). This sets the forward messages (`marginal_outgroup_LH`) of every node
2. `optimize_gtr_rate()` runs only the postorder pass in its cost function [packages/legacy/treetime/treetime/treeanc.py#L1685-L1689](../../packages/legacy/treetime/treetime/treeanc.py#L1685-L1689), so it updates the backward messages (`marginal_subtree_LH`) and leaves the forward messages unchanged
3. `infer_gtr()` skips the reconstruction when a marginal reconstruction already exists [packages/legacy/treetime/treetime/treeanc.py#L1544](../../packages/legacy/treetime/treetime/treeanc.py#L1544), and counts with both message directions [packages/legacy/treetime/treetime/treeanc.py#L1144-L1145](../../packages/legacy/treetime/treetime/treeanc.py#L1144-L1145)

v0 therefore also keeps its forward messages fixed during the later iterations, but they come from the first fitted GTR, not from the initial uniform model.

## Impact

The divergence affects every iteration after the first, and `--iterations` defaults to 5 [packages/app-commands/src/commands/mugration/args.rs#L143-L149](../../packages/app-commands/src/commands/mugration/args.rs#L143-L149), so default runs are affected. The size of the effect is not measured. It is a candidate cause of the confidence differences from v0 in [H-mugration-iterative-gtr.md](H-mugration-iterative-gtr.md), which reach 3.5e-2 on `dengue_20` and change trait assignments at near-tied nodes.

> [!IMPORTANT]
> **Investigation required.** Measure how much of the golden-master difference of [H-mugration-iterative-gtr.md](H-mugration-iterative-gtr.md) this divergence explains: implement the v0 order and rerun the ignored tests of [packages/treetime/src/mugration/__tests__/test_gm_mugration.rs](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs).

## Fix

Match v0: after the first GTR fit, run one full marginal update under the fitted model, and count the later iterations with its forward messages. This follows the reference, so it needs no decision.

[kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md) describes a further change, a full marginal update in every iteration (expectation maximization with consistent messages). That change diverges from v0 and needs a decision. Its section "Current behavior" states that v1 matches v0 here, which this issue shows is not the case.

The ancestral refinement of `--gtr-iterations` (`fn refine_gtr_model()`) is a separate case: v0's `infer_gtr_iterative()` runs a full reconstruction in each iteration, and [M-ancestral-gtr-iterations-refit-against-frozen-messages.md](M-ancestral-gtr-iterations-refit-against-frozen-messages.md) covers that divergence.

## Validation

- A test on a small tree compares the forward messages that the second iteration counts with against a marginal update under the first fitted GTR
- The ignored mugration golden-master tests are rerun, and their remaining difference is recorded in [H-mugration-iterative-gtr.md](H-mugration-iterative-gtr.md)

## Related issues

- [H-mugration-iterative-gtr.md](H-mugration-iterative-gtr.md)
- [M-mugration-gtr-rate-restore-inconsistency.md](M-mugration-gtr-rate-restore-inconsistency.md)
- [M-ancestral-gtr-iterations-refit-against-frozen-messages.md](M-ancestral-gtr-iterations-refit-against-frozen-messages.md)
