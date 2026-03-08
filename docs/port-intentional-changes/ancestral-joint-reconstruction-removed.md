# Joint ML ancestral reconstruction removed

| Property    | Value                                                                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Intentional feature removal                                                                                                                                                   |
| v1 Location | `MethodAncestral::Joint` (`#MethodAncestral`) in [`packages/treetime/src/commands/ancestral/args.rs#L11-L16`](../../packages/treetime/src/commands/ancestral/args.rs#L11-L16) |
| v1 Runtime  | `unimplemented!()` at [`packages/treetime/src/commands/ancestral/run.rs#L193-L194`](../../packages/treetime/src/commands/ancestral/run.rs#L193-L194)                          |
| v0 Location | `_ml_anc_joint()` (`#_ml_anc_joint`) in [`packages/legacy/treetime/treetime/treeanc.py#L934-L1080`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1080)            |
| Affects     | `--method-anc joint` (v0 default and v1 default, panics at runtime in v1)                                                                                                     |
| Commands    | `ancestral`, `timetree`, `clock`                                                                                                                                              |

## What joint reconstruction does

Joint ML reconstruction (Pupko et al. 2000) finds the single most probable
global assignment of ancestral states across the entire tree. It uses a
Viterbi-style algorithm:

- Postorder pass: each node accumulates `joint_Lx` (log-likelihood of best
  subtree assignment) and `joint_Cx` (argmax backpointers) via GTR log-transition
  matrices.
- Root: selects the most probable state at each position using
  `joint_Lx + log(pi)`.
- Preorder pass: traces back through `joint_Cx` pointers to reconstruct the
  globally optimal sequence at each internal node.

v0 implementation
([`treeanc.py#L960-L1063`](../../packages/legacy/treetime/treetime/treeanc.py#L960-L1063))
operates entirely in log-space, accumulating via addition.

## Why v1 removes it

Marginal reconstruction computes per-site posterior probabilities independently
at each node. Joint reconstruction picks the single most likely global
assignment, collapsing probability distributions into point estimates. This
discards per-site uncertainty information that downstream steps (GTR inference,
timetree optimization, convergence monitoring) rely on.

v1's architecture is built around probability profiles (dense or sparse) that
flow through the tree during belief propagation. Joint reconstruction produces
sequences, not profiles, making it incompatible with the profile-based pipeline.

The `MethodAncestral::Joint` enum variant is preserved (and is `#[default]`) to
maintain CLI argument compatibility with v0 invocations. It panics at runtime
with `unimplemented!()`.

## Practical impact

- `treetime ancestral` without `--method-anc` panics (joint is default).
- Users must specify `--method-anc marginal` or `--method-anc parsimony`.
- No golden master tests exist for joint reconstruction in v1.
- Marginal reconstruction produces strictly more information (full posteriors)
  and is the recommended method for all use cases.
