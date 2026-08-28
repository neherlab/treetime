# Parallelize multi-partition marginal reconstruction

## Summary

`treetime ancestral` reconstructs each per-CDS amino-acid translation as an independent marginal partition on the shared tree. These partitions are now processed one at a time (a landed memory fix). Each partition is independent, so the outer loop over partitions is a candidate for parallel execution. This proposal decomposes the parallelization decision into orthogonal axes, records the trade-offs per axis, and recommends a combination to implement. It also generalizes the design so future multi-partition work for nucleotide sequences and scalar traits reuses the same driver.

## Context and background

- The multi-partition machinery lives in `packages/treetime/src/ancestral/multi.rs`. `fn reconstruct_marginal_partition()` [packages/treetime/src/ancestral/multi.rs#L56](packages/treetime/src/ancestral/multi.rs#L56) builds one partition, runs the marginal passes, reconstructs node states, and returns the result. `struct PartitionPlan` carries an `alphabet` field, so the function is alphabet-agnostic, not amino-acid specific.
- The only caller today is the amino-acid path: the per-CDS loop `for (index, cds) in cdses.iter().enumerate()` [packages/treetime/src/commands/ancestral/run.rs#L311](packages/treetime/src/commands/ancestral/run.rs#L311). Per CDS it:
  - reads one translation FASTA,
  - builds a `PartitionPlan` and reconstructs,
  - extracts `AaCdsNodeData` and writes the per-CDS FASTA,
  - drops the partition before the next CDS.
- Measured on the `nextstrain/mpox` `hmpxv1` bundle (179 CDSes): holding every partition resident used about 8.7 GB peak RSS; processing one at a time uses about 0.74 GB, byte-identical output. Each resident partition costs about 45 MB (per-edge probability vectors over the roughly 20-symbol amino-acid alphabet).

### Verified findings that shape the design

- F1 (inner parallelism already exists). The sparse marginal backward and forward passes run through an indexed pass whose workers process nodes as their child-to-parent dependencies clear, inside a `rayon::scope` worker pool [packages/treetime/src/partition/indexed_pass.rs#L437](packages/treetime/src/partition/indexed_pass.rs#L437). Worker count is `rayon::current_num_threads()` [packages/treetime/src/partition/indexed_pass.rs#L391](packages/treetime/src/partition/indexed_pass.rs#L391), so inner parallelism sizes to the current rayon pool. Each node's math is local and the traversal is dependency-ordered, so results are numerically deterministic regardless of worker count.
- F2 (RNG is confined to the sampler). Only `fn ancestral_reconstruction_marginal()` [packages/treetime/src/ancestral/marginal.rs#L43](packages/treetime/src/ancestral/marginal.rs#L43) consumes the random number generator, and only under `--sample-from-profile=root|all`. The expensive phase, `fn update_marginal()` [packages/treetime/src/ancestral/marginal.rs#L30](packages/treetime/src/ancestral/marginal.rs#L30) (backward and forward passes) and per-partition GTR inference, uses no RNG. The default `argmax` mode draws nothing and is fully deterministic.
- F3 (discrete-trait partitions do not fit the current driver bound). The driver is bound on `trait MarginalAugurPartition` [packages/treetime/src/ancestral/multi.rs#L117](packages/treetime/src/ancestral/multi.rs#L117), which requires `PartitionMarginalOps` and `AugurNodeDataJsonAncestralPartition`. The scalar-trait partition `struct PartitionMarginalDiscrete` implements `PartitionMarginalPasses` [packages/treetime/src/partition/marginal_discrete.rs#L204](packages/treetime/src/partition/marginal_discrete.rs#L204) but not the ancestral-sequence traits; mugration extracts state-confidence maps, not sampled sequences. `trait PartitionMarginalOps` extends `trait PartitionMarginalPasses` [packages/treetime/src/partition/traits.rs#L209](packages/treetime/src/partition/traits.rs#L209), so `PartitionMarginalPasses` [packages/treetime/src/partition/traits.rs#L197](packages/treetime/src/partition/traits.rs#L197) is the common denominator for the parallel phase.

## Goals

- Reduce wall-clock time of multi-partition reconstruction by running independent partitions in parallel.
- Preserve a deterministic mode that produces byte-identical output to the current sequential result.
- Keep peak memory bounded and never oversubscribe CPU cores.
- Make the driver generic so nucleotide multi-partitioning and scalar-trait partitions reuse it.

## Non-goals

- Changing the numerical algorithm or the argmax reconstruction output.
- Parallelizing the single nucleotide main reconstruction that currently runs through `pipeline::run` (separate path; benefits only after nucleotide multi-partitioning routes through this driver).
- Exact-reproducible parallel sampling (a counter-based-RNG option is recorded but not recommended now).

## Design axes

### Axis A: outer parallelism model

- A1. Sequential loop (current). One partition at a time.
  - Pros: simplest, lowest memory (about one partition resident), deterministic.
  - Cons: single core for the outer loop; no throughput gain across partitions.
- A2. **(Recommended)** Bounded outer concurrency with rayon over partitions. Process up to `OUTER` partitions concurrently.
  - Pros: near-linear throughput across partitions where inner parallelism is not already saturating cores; memory bounded by `OUTER`.
  - Cons: overlapping partitions make sampling draw order non-deterministic (see Axis C).
- A3. Unbounded `par_iter` over all partitions.
  - Pros: trivial to write.
  - Cons: resident partitions scale with concurrency up to the CDS count, reintroducing the memory blow-up the landed fix removed. Rejected.

### Axis B: thread-budget knobs

- B1. Single knob `-j` (total threads), outer concurrency fixed.
  - Pros: one number to reason about.
  - Cons: cannot separate partition concurrency (memory and determinism) from per-partition CPU.
- B2. **(Recommended)** Two knobs: `-j` = total thread budget; `--jobs-for-partitions` = `OUTER` (outer concurrency), with `1 <= OUTER <= -j`. Inner threads per partition is about `-j / OUTER`.
  - Pros: `OUTER` becomes an independent memory-and-determinism knob; `OUTER=1` yields deterministic output while inner threads still use all `-j` cores (F1); `-j` caps total CPU.
  - Cons: two numbers; requires the wiring rule in Axis E to avoid the inner pass reading the wrong thread count.
  - Oversubscription rule: clamp so `OUTER * inner <= -j <= cores`. Expressing both knobs as a factorization of the budget makes oversubscription structurally impossible.

### Axis C: RNG and determinism strategy

Applies only to `--sample-from-profile=root|all`; argmax draws nothing and is always deterministic (F2).

- C1. **(Recommended)** Per-partition RNG, non-deterministic under `OUTER > 1`. `OUTER=1` (or `-j 1`) processes partitions in fixed order and is byte-identical to the current sequential output.
  - Pros: simplest to build and maintain; deterministic mode available on demand; no change to the numerical core.
  - Cons: sampled sequences vary run to run when `OUTER > 1`. Requires documenting the contract and recording a decision.
- C2. Per-partition seeded RNG from `(seed, cds_name)`, deterministic at any `OUTER`.
  - Pros: reproducible even in parallel.
  - Cons: changes the exact sampled bytes versus today (independent streams instead of one ordered stream), so it is a one-time reproducibility-baseline change requiring consent. Needs a little more plumbing than C1.
- C3. Counter-based seekable RNG (for example ChaCha with word-position seek), exact-parallel. Precount draws per node, prefix-sum to per-draw offsets, sample in parallel by seeking.
  - Pros: full parallelism with output identical to a single ordered stream.
  - Cons: much more complex. Each logical draw must consume a fixed number of RNG words, which requires verifying the categorical sampler is inverse-CDF rather than rejection based. Matching today's exact bytes also requires the current RNG to be seekable and fixed-width. Rejected now as over-engineered for the benefit.

### Axis D: driver location and generality

- D1. Parallelize inside the amino-acid loop in `run.rs`.
  - Pros: smallest diff.
  - Cons: amino-acid only; future nucleotide-multi and scalar-trait partitions each re-implement parallelism and the thread and memory knobs.
- D2. **(Recommended)** Generic partition-parallel driver in `multi.rs` over a `Vec<PartitionPlan>` plus a pluggable per-partition consumer closure. Bind on `PartitionMarginalPasses` (the common denominator, F3) for the parallel phase; move output assembly (amino-acid node data, trait confidence maps, nucleotide mutations) into the consumer; keep the RNG-driven sequence sampling inside the consumer so trait partitions that never sample are unaffected.
  - Pros: one implementation serves amino-acid, future nucleotide-multi, and scalar-trait partitions; knobs and determinism contract defined once.
  - Cons: requires widening the trait bound and routing a consumer type through the driver.

### Axis E: inner parallelism handling (multiple choice)

- E1. **(Recommended)** Reuse the existing indexed-pass inner parallelism (F1) for the sparse path; add no new inner code. Wire each partition's prepare phase to run under the inner pool so `rayon::current_num_threads()` reports `INNER`, not `OUTER` (otherwise the inner pass spawns `OUTER` workers into the outer pool and steals the outer budget).
- E2. **(Recommended, conditional)** Confirm and, if absent, add equivalent inner parallelism to the dense and discrete indexed passes, so scalar-trait and nucleotide-dense partitions also benefit at `OUTER=1`. Keep any inner reductions order-stable to preserve numerical determinism.

## Recommended combination

The recommended combination across axes:

- A2 bounded outer concurrency.
- B2 two knobs: `-j` (total budget) and `--jobs-for-partitions` (`OUTER`).
- C1 per-partition non-deterministic RNG with a deterministic `OUTER=1` mode.
- D2 generic driver in `multi.rs`.
- E1 reuse existing inner parallelism with the `inner_pool.install` wiring rule.

What this achieves:

- Partitions run in parallel with per-partition CPU still parallelized internally.
- `OUTER` bounds both memory (about `OUTER * 45` MB) and determinism.
- Oversubscription is impossible by the clamp.
- The driver serves all partition types.
- The only behavioral change is non-deterministic sampled sequences under `OUTER > 1`, which is opt-in and documented.

Defer C3 (exact-parallel sampling) and treat C2 as the fallback only if a decision requires deterministic parallel sampling.

## Determinism contract

| Mode                                  | OUTER = 1                    | OUTER > 1                           |
| ------------------------------------- | ---------------------------- | ----------------------------------- |
| argmax (default)                      | deterministic                | deterministic                       |
| `--sample-from-profile=root` or `all` | byte-identical to sequential | non-deterministic sampled sequences |

- Inner thread count never affects results (F1): the passes are dependency-ordered with node-local math.
- Document in the `--sample-from-profile` help text: reproducible sampling requires `--jobs-for-partitions=1` (or `-j 1`); argmax is unaffected.

## Memory and thread model

- Resident partitions equal `OUTER`; peak footprint is about `OUTER * per-partition-size` (about 45 MB for amino-acid CDSes; larger for long nucleotide partitions; small for scalar traits).
- Total worker threads equal `OUTER * INNER`; clamp `OUTER * INNER <= -j <= cores`.
- `OUTER=1` gives the minimum memory and the deterministic mode while still using `-j` cores through inner parallelism.

## Generalization to other partition types

- Nucleotide multi-partitioning (future): route through the same driver by supplying nucleotide `PartitionPlan`s and a nucleotide consumer. The current single nucleotide reconstruction through `pipeline::run` is separate and out of scope here.
- Scalar traits (future multi-trait mugration): `PartitionMarginalDiscrete` already satisfies `PartitionMarginalPasses` (F3); it fits the driver once the bound is lowered and output assembly moves to the consumer. Trait partitions never sample sequences, so the RNG contract does not apply to them.

## Open questions to resolve before or during implementation

- Q1. Confirm whether the dense and discrete indexed passes are actually parallel (only the sparse path is verified). This decides whether scalar-trait and nucleotide-dense partitions benefit from `OUTER=1` or need inner parallelism added (Axis E2).
- Q2. If C2 or C3 is ever chosen, verify the categorical sampler's RNG word consumption is fixed per draw (inverse-CDF) rather than rejection based.
- Q3. Decide default values: proposed `-j` defaults to the core count and `--jobs-for-partitions` defaults to 1 (deterministic and lowest memory by default; users opt into `OUTER > 1` for throughput). Confirm this default policy.

## Validation plan and acceptance criteria

- Equivalence: `-j 1` output byte-identical to the pre-change output for both argmax and sampling.
- Parallel argmax: `-j N` with any `OUTER` byte-identical to `-j 1` (deterministic).
- Deterministic sampling mode: `--jobs-for-partitions=1` with `-j N` byte-identical to sequential; two runs identical.
- Non-deterministic sampling: `OUTER > 1` produces valid output (structure and per-node fields well-formed); reproducibility not asserted.
- No oversubscription: assert `OUTER * INNER <= -j` for all knob combinations.
- Memory: peak RSS scales with `OUTER`, not with partition count.
- Set pool sizes per test via a scoped rayon `ThreadPool`.

## Follow-up

- Extract the recommended combination into `kb/issues/` as separate actionable items (driver refactor and bound widening; two-knob CLI and clamp; per-partition RNG and determinism contract; documentation; the dense/discrete inner-parallelism verification from Q1). Create tickets only for issues whose design axes above are decided.
- On implementation, record the non-deterministic-sampling divergence in `kb/decisions/` (requires explicit human consent, since it changes sampling reproducibility under parallel execution).
