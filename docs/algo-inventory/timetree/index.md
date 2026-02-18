# Timetree Inference Algorithms

[Back to index](../index.md)

## Belief Propagation

| Property    | Value                                                                                   |
| ----------- | --------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                              |
| v1 Location | `packages/treetime/src/commands/timetree/inference/backward_pass.rs`, `forward_pass.rs` |
| v0 Location | `packages/legacy/treetime/treetime/node_interpolator.py`                                |
| Reference   | Pearl, J. (1988). "Probabilistic Reasoning in Intelligent Systems." Morgan Kaufmann     |

Two-pass message passing for time inference: backward (convolve + multiply child distributions), forward (divide + convolve parent distributions).

---

## Kingman Coalescent

| Property    | Value                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                     |
| v1 Location | `packages/treetime/src/commands/timetree/coalescent/*.rs` (6 files)                            |
| v0 Location | `packages/legacy/treetime/treetime/merger_models.py`                                           |
| Reference   | Kingman, J.F.C. (1982). "The coalescent." Stochastic Processes and Applications, 13(3):235-248 |

Computes lineage dynamics k(t), integral merger rate I(t), per-node survival/merger contributions.

---

## Tc Optimization

| Property    | Value                                                                   |
| ----------- | ----------------------------------------------------------------------- |
| Type        | Well-known optimization                                                 |
| v1 Location | `packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs:43-` |
| Reference   | Brent 1973                                                              |

Optimizes coalescent time scale in log space [-20, 2].

---

## Skyline Coalescent

| Property    | Value                                                                                         |
| ----------- | --------------------------------------------------------------------------------------------- |
| Type        | Well-known (Nelder-Mead) + domain-specific (skyline)                                          |
| v1 Location | `packages/treetime/src/commands/timetree/coalescent/skyline.rs:68-`                           |
| v0 Location | `packages/legacy/treetime/treetime/merger_models.py:281-`                                     |
| Reference   | Strimmer & Pybus (2001). "Exploring the demographic history." Mol Biol Evol, 18(12):2298-2305 |
| Reference   | Minin et al. (2008). "Smooth skyride through a rough skyline." Mol Biol Evol, 25(7):1459-1471 |

Piecewise-varying Tc estimation with smoothness regularization.

**Modification**: v1 uses Nelder-Mead; v0 uses SLSQP.

---

## Relaxed Clock

| Property    | Value                                                                                                                        |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                   |
| v1 Location | `packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs:25-`                                                  |
| Reference   | Thorne et al. (1998). "Estimating the rate of evolution of the rate of molecular evolution." Mol Biol Evol, 15(12):1647-1657 |

Two-pass quadratic penalty optimization for rate variation.

---

## Polytomy Resolution

| Property    | Value                                                                  |
| ----------- | ---------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                             |
| v1 Location | `packages/treetime/src/commands/timetree/optimization/polytomy.rs:27-` |

Greedy pairwise merging with Brent optimization per pair.

**Note**: Stochastic resolution (v0's `generate_subtree()`) not ported. See [unimplemented](../unimplemented/#stochastic-polytomy-resolution).

---

## Iterative EM-like Refinement

| Property    | Value                                                                          |
| ----------- | ------------------------------------------------------------------------------ |
| Type        | Well-known (EM)                                                                |
| v1 Location | `packages/treetime/src/commands/timetree/refinement.rs:21-`                    |
| Reference   | Sagulenko et al. (2018). "TreeTime." Virus Evolution, 4(1):vex042, Section 2.4 |

Alternates sequence reconstruction (E-step) and time inference (M-step).

---

## File Index

| File                                                        | Algorithms                                   |
| ----------------------------------------------------------- | -------------------------------------------- |
| `packages/treetime/src/commands/timetree/inference/*.rs`    | Belief propagation, branch distributions     |
| `packages/treetime/src/commands/timetree/coalescent/*.rs`   | Kingman coalescent, skyline, Tc optimization |
| `packages/treetime/src/commands/timetree/optimization/*.rs` | Polytomy, relaxed clock, reroot              |
