# Command relationships: prune, optimize, timetree

Scientific and architectural relationships between the three main tree-refinement commands. Covers current implementation, ideal design, the two EM-like loops, and a gap table mapping known issues to required changes.

## Current implementation: Venn diagram

```mermaid
graph TD
    subgraph prune_only["prune only"]
        P1[topology collapse]
        P2[shared-mutation polytomy merge]
        P3[prune by name / length / empty]
    end

    subgraph prune_optimize_shared["prune ∩ optimize"]
        PO1[Fitch compression]
        PO2[sparse partition]
    end

    subgraph optimize_only["optimize only"]
        O1[Newton/grid branch length optimization]
        O2[exponential damping]
        O3[zero-branch derivative-sign detection]
    end

    subgraph optimize_timetree_shared["optimize ∩ timetree"]
        OT1[GTR inference]
        OT2[marginal reconstruction — update_marginal]
        OT3[EM-like iterative loop]
    end

    subgraph timetree_only["timetree only"]
        T1[clock regression / reroot]
        T2[node time belief propagation]
        T3[date constraints]
        T4[coalescent prior / skyline]
        T5[confidence intervals]
        T6[relaxed clock]
        T7[temporal polytomy resolution]
    end

    subgraph all_shared["all three"]
        A1[tree input]
        A2[fasta input]
        A3[alphabet]
    end

    prune_only --- PO1
    optimize_only --- PO1
    optimize_only --- OT1
    timetree_only --- OT1
```

## Current key intersections

| Intersection | Shared scientific components |
|:-------------|:-----------------------------|
| prune ∩ optimize | Fitch-compressed sparse partition (structure, no inference) |
| optimize ∩ timetree | GTR inference + `update_marginal` ancestral reconstruction (shared code, different M-step objective) |
| prune ∩ timetree | — (none beyond I/O) |
| all three | tree input, fasta input, alphabet |

## The two EM-like loops compared

Both optimize and timetree use an EM-like alternating optimization. The E-step is identical. The M-steps differ fundamentally.

| Phase | optimize | timetree |
|:------|:---------|:---------|
| E-step | `update_marginal()` — ancestral reconstruction | `update_marginal()` — same call, same code |
| M-step | `run_optimize_mixed()` — Newton/grid on **branch lengths** (free, unconstrained) | `run_timetree()` — belief propagation on **node times** (constrained by clock, dates, coalescent) |
| Primary variable | branch length $b_e$ (substitution units) | node time $t_v$ (calendar units) |
| Constraint | none — free $b_e \geq 0$ | clock: $b_e = \Delta t \cdot \mu \cdot \gamma$; date bounds; coalescent prior |
| Topology cleanup | absent (known issue `M-optimize-no-topology-cleanup-in-loop`) | temporal greedy polytomy resolution |
| Damping | explicit post-sweep blending (`d = 0.75`) | absent |
| Convergence | $\|\Delta \ell\| < \epsilon$ | `n_diff == 0 && n_resolved == 0` |

The two loops are **siblings** sharing an E-step, not parent/child. Running optimize's M-step inside timetree's loop would optimize the wrong objective — unconstrained per-edge substitution rate when the constraint is a molecular clock.

## Current objective functions

| Command | Primary optimization variable | Objective |
|:--------|:------------------------------|:----------|
| prune | — | topology manipulation only |
| optimize | branch lengths (substitution space) | $\ell(\text{branch lengths} \mid \text{sequences, GTR})$ |
| timetree | node times (calendar space) | $\ell(\text{node times} \mid \text{sequences, GTR, clock rate, dates, coalescent})$ |

---

## Ideal implementation

### The ideal relationship hierarchy

```
topology_cleanup/               ← shared module, no GTR, no ML
  collapse_sparse_edge()        ← composition-correct (not union)
  prune_short_branches()
  merge_shared_mutation_branches()

prune/
  run_prune()
    topology_cleanup::prune_nodes()      ← by name / length / empty
    topology_cleanup::merge_shared()     ← optional

optimize/
  run_optimize()
    initial_guess_mixed()
    loop:
      update_marginal()                  ← E-step
      run_optimize_mixed()               ← M-step
      apply_damping(0.75)
      topology_cleanup::prune_short()    ← NEW: inside loop, after damping
      topology_cleanup::merge_shared()   ← NEW: after prune exposes polytomies
      graph.build()                      ← only if topology changed

timetree/
  run_timetree_estimation()
    optimize(max_iter=1)                 ← NEW: pre-step, seeds branch lengths
    clock_regression + reroot
    loop:
      update_marginal()                  ← E-step (shared with optimize)
      resolve_polytomies_temporal()      ← M-step topology (time-domain, distinct)
      run_timetree()                     ← M-step time inference
      coalescent / relaxed_clock
```

### Why prune ⊆ optimize (topology cleanup inside the loop)

EM theory (Dempster et al. 1977) guarantees monotonic likelihood increase only when the E-step and M-step operate on a consistent topology. Zero-length edges left in the tree cause the E-step to reconstruct ancestors for nodes that should not exist, and the M-step to optimize branch lengths for phantom edges. **Pruning inside the loop is required for the EM convergence argument to hold**, not just a quality-of-life feature.

Order within topology cleanup: **prune first, then merge**. Merging before pruning hides shared mutations behind short internal nodes that pruning would have exposed as part of the same polytomy (`M-prune-wrong-operation-order`).

The `prune` command as a user-facing entry point retains value — users want topology cleanup without full ML. Its implementation should delegate to the same shared `topology_cleanup` module used inside optimize's loop (Chapter 10 of the iterative tree refinement report, Recommendation R3).

### Why optimize(max_iter=1) precedes timetree's loop

Time inference uses branch-length distributions $P(b_e) = \exp(Q b_e)$ to inform node time posteriors. If the input tree carries raw neighbor-joining or parsimony branch lengths, these distributions are poorly centered so the initial time posteriors are inaccurate. One ML optimization pass seeds the loop with branch lengths already close to the ML optimum (`M-timetree-missing-initial-branch-optimization`). This matches v0's design exactly (`treetime.py:243, 266`).

---

## Current vs ideal: gap table

| Gap | Severity | Tracking issue | Blocker |
|:----|:---------|:---------------|:--------|
| `collapse_sparse_edge` uses `iterator_union` instead of `compose_substitutions` | Blocking | `M-prune-collapse-uses-union-not-composition` | Prerequisite for all topology integration |
| Optimize loop has no topology cleanup (zero-length edges accumulate) | High | `M-optimize-no-topology-cleanup-in-loop` | Needs composition fix first |
| Prune runs merge before prune (wrong order) | Medium | `M-prune-wrong-operation-order` | Standalone fix |
| Topology cleanup not in a shared module | Medium | Chapter 10 (iterative-tree-refinement), R3 | None |
| Timetree skips optimize(max_iter=1) pre-step | Medium | `M-timetree-missing-initial-branch-optimization` | None |
| Timetree has no branch re-optimization after polytomy resolution | Low | Chapter 10 (iterative-tree-refinement) | None |

---

## Ideal Venn diagram

In the ideal design, prune is a proper subset of optimize (shared topology_cleanup module), and optimize provides a one-pass initialization to timetree before its own loop begins. The two loops remain siblings sharing an E-step.

```mermaid
graph TD
    subgraph topology_cleanup["topology_cleanup — shared module"]
        TC1[compose_substitutions edge collapse]
        TC2[prune_short_branches]
        TC3[merge_shared_mutation_branches]
    end

    subgraph prune_cmd["prune command"]
        P1[prune by name / length / empty]
    end

    subgraph optimize_loop["optimize loop"]
        O1[Newton/grid branch length M-step]
        O2[exponential damping]
        O3[zero-branch detection]
    end

    subgraph shared_estep["shared E-step"]
        E1[update_marginal — ancestral reconstruction]
        E2[GTR inference]
    end

    subgraph timetree_loop["timetree loop"]
        T1[temporal belief propagation M-step]
        T2[clock / date / coalescent constraints]
        T3[temporal polytomy resolution]
        T4[confidence intervals]
    end

    prune_cmd --> TC1
    prune_cmd --> TC2
    prune_cmd --> TC3
    optimize_loop --> TC1
    optimize_loop --> TC2
    optimize_loop --> TC3
    optimize_loop --> E1
    optimize_loop --> E2
    timetree_loop --> E1
    timetree_loop --> E2
    optimize_loop -- "pre-step max_iter=1" --> timetree_loop
```
