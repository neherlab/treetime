# Proposal: stochastic coalescent polytomy resolution in the timetree loop

## Summary

Replace the greedy pairwise merger in
[packages/treetime/src/timetree/optimization/polytomy/mod.rs](../../packages/treetime/src/timetree/optimization/polytomy/mod.rs)
with a stochastic coalescent-with-mutations sweep, modelled on v0's `generate_subtree`
(reached via `--stochastic-resolve`). The sweep runs backwards in time from the most recent
child of a polytomy toward the parent, drawing either a mutation event or a coalescence
event at each step. Only branches that are already alive and have placed all their mutations
may coalesce.

Three defects in the v0 loop are not reproduced; each has its own erratum entry. The greedy
implementation is removed rather than kept behind a flag.

## Current state

`fn resolve_polytomies()` scans for nodes with more than two children and merges pairs
greedily by likelihood gain
([packages/treetime/src/timetree/optimization/polytomy/mod.rs#L166-L318](../../packages/treetime/src/timetree/optimization/polytomy/mod.rs#L166-L318)).
Children are first partitioned into "stretched" (`mutation_length < clock_length`) and
"compressed"; only stretched children are merged by default. Each candidate pair is scored by
Brent-optimizing the merge time, and the highest-gain pair is merged until the gain falls
below `DEFAULT_RESOLUTION_THRESHOLD`.

This is v0's `_poly` path
([packages/legacy/treetime/treetime/treetime.py#L713-L870](../../packages/legacy/treetime/treetime/treetime.py#L713-L870)),
which v0 itself deprecates: `resolve_polytomies` emits a deprecation warning stating the
greedy mode "is not well suited for large polytomies" and that stochastic resolution "will
become the default in future versions"
([packages/legacy/treetime/treetime/treetime.py#L679-L690](../../packages/legacy/treetime/treetime/treetime.py#L679-L690)).

## Motivation

Greedy pairwise merging by maximum likelihood gain produces atypical topologies. It commits
to the single highest-scoring merger at each step, so the resulting subtree is a mode of the
scoring function rather than a draw from any process. On large polytomies this concentrates
structure into a few deep cherries and leaves the remainder unresolved when the gain drops
below threshold.

The stochastic sweep instead samples from a coalescent process conditioned on the observed
mutations: a child carrying $m$ substitutions cannot coalesce until all $m$ have been placed
in time, so mutation-laden branches stay separate longer and mutation-free branches merge
early. The resulting subtree is a plausible realisation rather than an extremum, and the
number of children actually resolved falls out of the available time window rather than a
likelihood threshold.

## Algorithm

At a polytomy with parent $P$ and children $c_1 \dots c_k$:

- $m_b$ -- substitutions remaining to be placed on branch $b$, initially the count on the
  child edge
- `alive` -- branches the sweep has reached; `to_come` -- those it has not
- `ready` $= \{b \in \text{alive} : m_b = 0\}$ -- the only branches eligible to coalesce

The sweep starts at the most recent child and moves toward the parent. At each step it draws
a waiting time from the total event rate, then either removes one mutation from a branch or
merges two `ready` branches under a new node. It terminates when two lineages remain or the
parent's time is reached; survivors attach to the parent as a residual polytomy.

### Time convention

This is the main translation hazard. v0 works in `time_before_present`, increasing into the
past; v1 uses calendar time with `parent.time < child.time`. Every comparison inverts:

| v0                                               | v1                                             |
| ------------------------------------------------ | ---------------------------------------------- |
| `tmax = parent.time_before_present`              | `t_stop = parent.time` (lower bound)           |
| sort children ascending by `time_before_present` | sort **descending** by `time`                  |
| `t` starts at most recent child, increases       | `t` starts at `max(child.time)`, **decreases** |
| `while t < tmax`                                 | `while t > t_stop`                             |
| `t += dt`                                        | `t -= dt`                                      |
| pop when `t > to_come[0].tbp`                    | pop when `t < to_come[0].time`                 |
| `b.branch_length = tmax - b.time`                | `time_length = b.time - t_stop`                |
| early return `if t >= tmax`                      | early return if `max(child.time) <= t_stop`    |

### Rates

Write $\mu$ for the whole-alignment mutation rate, $\kappa(t)$ for the per-branch coalescent
merger rate, $M = \sum_b m_b$, $R$ for `ready`. The adopted rates are

$$R_{\text{mut}} = \mu M, \qquad R_{\text{coal}} = \max(0, \lvert R\rvert - 1)\,\kappa(t)$$

Draw $\Delta t \sim \operatorname{Exp}(R_{\text{mut}} + R_{\text{coal}})$; take a mutation
event with probability $R_{\text{mut}}/(R_{\text{mut}} + R_{\text{coal}})$; pick the mutating
branch $\propto m_b$; pick the coalescing pair uniformly from $R$.

These differ from v0, which adds $\mu$ into both channels and produces branch-selection
weights inconsistent with its own mutation rate. See
[kb/v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md](../v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md).

Both rates are chosen by the caller and handed to `resolve_polytomies`, which holds no policy
about where either comes from: $\mu$ as a scalar, $\kappa$ as a function of calendar time.

$\mu$ is already computed in the caller: `clock_model.clock_rate() * total_length` is passed
as `zero_branch_slope` at
[packages/treetime/src/timetree/refinement.rs#L96](../../packages/treetime/src/timetree/refinement.rs#L96).
Rename the parameter `mutation_rate`; no new plumbing.

$\kappa(t)$ is `CoalescentModel::branch_merger_rate` on the same model the round's node times are
inferred under: lineage counts frozen before the optimization loop, $T_c$ as re-estimated for that
round. Assembled in
[packages/treetime/src/timetree/pipeline.rs](../../packages/treetime/src/timetree/pipeline.rs);
see
[kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md).
The frozen counts mean the sampler is not chasing a rate that drifts with the trees its own edits
produce, while $T_c$ still adapts. When the run carries no coalescent prior of its own the model
still exists, built from a constant $T_c$ estimated by the same one-segment analytic solve
`--coalescent-opt` uses.

v0 instead falls back to a dummy rate $\kappa = \tfrac12 \lvert R\rvert c$ with
$c = 2/(t_{\text{start}} - t_{\text{stop}})$ fixed at sweep start
([packages/legacy/treetime/treetime/treetime.py#L896](../../packages/legacy/treetime/treetime/treetime.py#L896)),
calibrated so the lineages typically coalesce within the window available above that one
polytomy. v1 does not reproduce it: it depends on both the polytomy's own window and the live
ready count, so it is not a function of time, and an estimated $T_c$ is a better-motivated
timescale than one manufactured from the window the answer has to fit into.

### Degenerate steps

When $R_{\text{mut}} + R_{\text{coal}} = 0$ -- no mutations among the alive branches and fewer
than two ready to coalesce -- no event can occur. If `to_come` is non-empty, advance `t` to the
next arrival time and re-enter. If it is empty, terminate; v0 relies on `Exp(0)` returning
infinity and the loop condition catching it
([packages/legacy/treetime/treetime/treetime.py#L922-L926](../../packages/legacy/treetime/treetime/treetime.py#L922-L926)),
which v1 should make explicit.

### Defects not reproduced

- Events committed past the parent bound, yielding negative branch lengths --
  [kb/v0-errata/timetree-stochastic-resolve-event-past-parent.md](../v0-errata/timetree-stochastic-resolve-event-past-parent.md).
  v1 tests the candidate time against `t_stop` before committing.
- The interval between a branch arrival and an overshooting draw is skipped --
  [kb/v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md](../v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md).
  v1 resumes at the arrival time.
- The rate/selection mismatch above.

## Bounded at the parent

The sweep hard-stops at `parent.time`; surviving lineages remain a residual polytomy. Pushing
further back would require relocating the parent node itself and re-checking against the
grandparent branch, i.e. mutating a second node's time inside a routine scoped to one node.

Deliberately out of scope: `run_timetree` re-infers every node time immediately afterwards
([packages/treetime/src/timetree/refinement.rs#L110](../../packages/treetime/src/timetree/refinement.rs#L110)),
so the parent's time is not final at this point. A partially resolved polytomy gets another
attempt on the next refinement iteration with a relocated parent, which reaches the same
outcome without coupling two nodes' times inside the sweep.

## Design

### Module layout

`timetree/optimization/polytomy/` (a directory, keeping `refinement.rs`'s import path valid):

- `mod.rs` -- public `resolve_polytomies`, plus the retained
  `validate_tree_before_topology_change`, `prepare_tree_after_topology_change`,
  `remove_single_child_nodes`, `inferred_time`
- `sweep.rs` -- the pure simulation
- `apply.rs` -- graph surgery

### Pure simulation

The sweep must not touch the graph. It takes a summary and returns a plan:

```rust
struct Lineage { id: usize, time: f64, mutations: u32 }
struct Merger  { time: f64, left: usize, right: usize }   // ids; new nodes get fresh ids
struct SubtreePlan { mergers: Vec<Merger>, roots: Vec<usize> }

fn simulate_subtree(
  children: &[Lineage], t_stop: f64, rates: &SweepRates, rng: &mut impl Rng,
) -> SubtreePlan;
```

`apply_plan(graph, parent_key, &plan)` then performs all surgery in one pass.

The split matters because the sweep is the part that is easy to get wrong and hard to debug
through a graph. As a pure function it is testable with a seeded RNG, property-testable, and
statistically testable against the Kingman waiting-time distribution without constructing a
tree at all.

### Mutation counts

Use `PartitionTimetreeOps::edge_subs(graph, edge_key)?.len()` summed over partitions. This is
exact, and better than v0's `round(mutation_length * L)`
([packages/legacy/treetime/treetime/treetime.py#L885](../../packages/legacy/treetime/treetime/treetime.py#L885)).
`subs_ml` is repopulated by the marginal forward pass
([packages/treetime/src/partition/marginal_passes.rs#L150](../../packages/treetime/src/partition/marginal_passes.rs#L150))
and `edge_subs()` errors when it is absent
([packages/treetime/src/partition/marginal_sparse.rs#L432](../../packages/treetime/src/partition/marginal_sparse.rs#L432)),
so fall back to the rounded product on error. Substitutions only, matching v0's
`mutation_length`; indels do not participate.

### Cost

The sweep runs $O(M + k)$ events. With an incrementally maintained `total_mutations` and a
swap-remove `Vec` for `ready`, per-event work is $O(1)$ except the weighted branch pick at
$O(\lvert W\rvert)$, giving $O(M \lvert W\rvert + k)$. A Fenwick tree over $m_b$ reduces the
pick to $O(\log k)$; deferred until profiling asks, but noted because large polytomies are the
stated motivation for the change.

### Supporting changes

**`Graph::reparent_edge(edge_key, new_source)`** in `treetime-graph` (~25 lines: pull the key
from the old source's `outbound_mut()`, push to the new, `set_source`). The existing
`fn merge_children()` reparents via `remove_edge` + `add_edge` with a fresh
`EdgeTimetree { time_length, ..default }`
([packages/treetime/src/timetree/optimization/polytomy/mod.rs#L539-L555](../../packages/treetime/src/timetree/optimization/polytomy/mod.rs#L539-L555)),
which discards the edge payload and changes the edge key. The graph operation must preserve the
edge key and payload because the sweep can reparent most children of a polytomy.

New edges get `time_length = child_time - parent_time`, `branch_length = Some(0.0)`,
`gamma = 1.0`.

**`CoalescentModel::branch_merger_rate(t)`** -- `compute_merger_rate_per_lineage_scalar`
([packages/treetime/src/coalescent/integration.rs#L12](../../packages/treetime/src/coalescent/integration.rs#L12))
is `pub(super)`. Expose it as a method mirroring the existing private `total_merger_rate`,
with the same finiteness and positivity guards
([packages/treetime/src/coalescent/coalescent.rs#L83-L100](../../packages/treetime/src/coalescent/coalescent.rs#L83-L100)).
Build the model once per refinement pass from `coalescent_tc` via `compute_coalescent_model`.

**RNG threading.** `Refinement` is constructed fresh inside the iteration loop
([packages/treetime/src/timetree/pipeline.rs#L309](../../packages/treetime/src/timetree/pipeline.rs#L309)),
so the generator must be created before the loop from `params.seed` via
`treetime_utils::sync::random::get_random_number_generator` and passed as
`&'a mut dyn RngCore` on the struct. That matches v0's single `self.rng` and keeps one
continuous stream across iterations.

`--seed` reaches `TimetreeParams` from
[packages/treetime/src/commands/timetree/args.rs](../../packages/treetime/src/commands/timetree/args.rs).
When it is `None`, the pipeline generates a seed, uses it, and logs it at `info` so the run can
be reproduced.

**Signature.** `resolve_polytomies` loses `resolution_threshold`, `merge_compressed` and
`clock_rate`; `zero_branch_slope` is renamed `mutation_rate`. `keep_polytomies` reaches
`TimetreeParams`
([packages/treetime/src/timetree/pipeline.rs#L61](../../packages/treetime/src/timetree/pipeline.rs#L61))
and is then unused, as tracked in
[kb/issues/N-timetree-unused-cli-flags.md](../issues/N-timetree-unused-cli-flags.md).

### Strategy selection

v1 has one stochastic path. v0 keeps a greedy path under `--greedy-resolve`, but v0
also deprecates it; carrying two resolution strategies would mean maintaining the scoring
machinery and its eight tests for a mode that is documented as unsuitable for the case it is
most often applied to.

## Consequences

- **Output is stochastic.** Runs without an explicit `--seed` are irreproducible; the
  logged seed is the recovery path.
- **No event-level v0 comparison.** Different RNG streams and corrected rates mean validation
  against v0 can only be distributional.
- **Convergence.** Once a polytomy is resolved into a binary subtree, later iterations find no
  polytomy there and do not re-randomise it, so the refinement loop still settles. Worth
  confirming empirically: resolution interacts with the convergence metrics recorded at
  [packages/treetime/src/timetree/pipeline.rs#L321-L329](../../packages/treetime/src/timetree/pipeline.rs#L321-L329).

## Validation plan

Pure sweep (`polytomy/__tests__/test_sweep.rs`), no graph:

- fixed seed produces a fixed plan (regression lock)
- invariants: every child appears exactly once as a plan leaf; merger times lie strictly
  inside $(t_{\text{stop}}, \min(\text{times of the merged pair}))$; a lineage never coalesces
  before it is alive; a lineage never coalesces with $m_b > 0$
- $\mu = 0$ with constant $\kappa$, many replicates: merger times match the Kingman
  waiting-time distribution (KS test)
- `max(child.time) <= t_stop` yields an empty plan
- zero total rate with `to_come` non-empty advances to the next arrival; with `to_come` empty
  it terminates rather than spinning
- the three errata: no merger is placed at $t < t_{\text{stop}}$; arrival resumption lands on
  the arrival time; the mutating-branch draw never falls through

Surgery and integration:

- a reparented child retains `branch_length` and `gamma`, and its sparse partition entry
  survives
- new edges carry `branch_length = 0.0`; all `time_length` values are positive
- no single-child nodes remain
- end-to-end refinement on a tree with a large polytomy: a seeded run is reproducible,
  converges, and leaves a sane node count

Datasets: ebola/20 and sc2/4500. Compare against v0 `--stochastic-resolve` distributionally
over replicates -- resulting child-count distribution per polytomy, tree-shape summaries,
final likelihood -- and against the current greedy implementation on the same trees to
characterise the change in resolved topology.

## Related

- [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md)
- [kb/v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md](../v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md)
- [kb/v0-errata/timetree-stochastic-resolve-event-past-parent.md](../v0-errata/timetree-stochastic-resolve-event-past-parent.md)
- [kb/v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md](../v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md)
- [kb/issues/N-timetree-unused-cli-flags.md](../issues/N-timetree-unused-cli-flags.md) --
  `--keep-polytomies`
- [kb/features/timetree.md](../features/timetree.md) -- polytomy resolution checklist
- [kb/decisions/optimize-polytomy-reversion-resolution.md](../decisions/optimize-polytomy-reversion-resolution.md) --
  the `optimize` loop's polytomy handling, which also uses `Graph::reparent_edge`
