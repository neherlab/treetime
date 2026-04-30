# Dense `edge_subs` optimization opportunities

Date: 2026-04-30
Base commit: `6396f86a` (`perf: change order of checks in dense edge_subs to run cheap checks first`)

## 1. Problem statement

Dense `edge_subs()` iterates every alignment position for each edge, calling `argmax_first()` on both parent and child profile rows to determine MAP states, then comparing them. Commit `6396f86a` reordered the checks to run cheap comparisons (`parent_state == child_state`, `is_canonical`) before `range_contains`.

This report measures one proposed follow-up (vectorized argmax) and catalogs further optimization opportunities for `edge_subs` and the surrounding dense code paths.

## 2. Setup

### 2.1 Datasets

Three datasets chosen to vary alignment length and tree size:

| Dataset           | Tips | Edges |  Alignment | Positions per call |
| ----------------- | ---: | ----: | ---------: | ------------------ |
| flu/h3n2/500      |  500 |  ~916 |   1,800 bp | 1.6M total         |
| mpox/clade-ii/20  |   20 |   ~39 | 197,000 bp | 7.7M total         |
| mpox/clade-ii/100 |  100 |  ~197 | 197,000 bp | 38.8M total        |

mpox datasets amplify alignment-length-sensitive costs (argmax, range_contains) because the alignment is 110x longer than flu.

### 2.2 Command

```bash
treetime ancestral --method-anc=marginal --dense=true \
  --tree=data/$v/tree.nwk --outdir=tmp/bench/ancestral \
  data/$v/aln.fasta.xz
```

The `ancestral` command calls `edge_subs()` once per edge during output annotation via `annotate_branch_mutations()`. Marginal reconstruction (forward/backward passes) dominates total runtime; `edge_subs` is a fraction.

### 2.3 Build

Release profile (LTO, opt-level 3) via Docker: `./dev/docker/run ./dev/dev br treetime`. Sequential OpenBLAS. Benchmarked with `hyperfine --warmup 1 --runs 5`.

### 2.4 Prior profile data

From the [dense OpenBLAS profiling report](dense-openblas-profiling.md) (sequential OpenBLAS, flu/h3n2/500):

| Category                       | Share | Functions                                                          |
| ------------------------------ | ----: | ------------------------------------------------------------------ |
| Math primitives (log/exp/norm) | 36.2% | `log_fma`, `softmax_with_log_norm`, `normalize_inplace`, `exp_fma` |
| GTR inference                  | 13.0% | `get_branch_mutation_matrix`, `accumulate_mutation_counts`         |
| Argmax / sequence extract      | 12.6% | `argmax_first` (all call sites, not just `edge_subs`)              |
| Other treetime                 |  3.3% | `assign_sequence`, `edge_subs`                                     |

`argmax_first` at 12.6% includes calls from marginal passes (`prof2seq`, `assign_sequence`, GTR inference) across all nodes. The `edge_subs` function itself falls within the 3.3% "other" category.

## 3. Experiment: vectorized argmax

Pre-compute argmax index vectors for entire parent and child profiles in a single pass, compare indices directly, and only perform character lookups and range checks at disagreeing positions.

### 3.1 Before (per-position scalar loop)

```rust
for (pos, parent, child) in izip!(0..nrows, parent_profile.rows(), child_profile.rows()) {
  let parent_state = alphabet.char(argmax_first(&parent).unwrap_or(0));
  let child_state = alphabet.char(argmax_first(&child).unwrap_or(0));
  if parent_state == child_state { continue; }
  if !alphabet.is_canonical(parent_state) || !alphabet.is_canonical(child_state) { continue; }
  if range_contains(parent_non_char, pos) || range_contains(child_non_char, pos) { continue; }
  subs.push(Sub::new(parent_state, pos, child_state)?);
}
```

Each iteration calls `argmax_first` twice (one per profile row), which internally loops over `n_states` elements (5 for nucleotides). The `ArrayView` row extraction from `rows()` involves per-row pointer arithmetic and bounds metadata.

### 3.2 After (vectorized argmax, index comparison)

```rust
let parent_indices = argmax_first_per_row(parent_profile);
let child_indices = argmax_first_per_row(child_profile);

for (pos, (&pi, &ci)) in parent_indices.iter().zip(child_indices.iter()).enumerate() {
  if pi == ci { continue; }
  let parent_state = alphabet.char(pi);
  let child_state = alphabet.char(ci);
  if !alphabet.is_canonical(parent_state) || !alphabet.is_canonical(child_state) { continue; }
  if range_contains(parent_non_char, pos) || range_contains(child_non_char, pos) { continue; }
  subs.push(Sub::new(parent_state, pos, child_state)?);
}
```

`argmax_first_per_row` iterates the `Array2` row-major memory in two contiguous passes (one per profile), producing `Vec<usize>` index vectors. The inner loop then compares `usize` pairs with no `ArrayView` overhead. Character lookups, canonical checks, and range checks only run at disagreeing positions (<5% on typical phylogenetic edges).

### 3.3 Test verification

2779 tests pass. 1 pre-existing known failure (`test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency`) is unrelated.

## 4. Results

| Dataset           | Baseline        | Optimized       | Wall speedup | User CPU saved |
| ----------------- | --------------- | --------------- | ------------ | -------------- |
| flu/h3n2/500      | 775 ms +/- 13   | 739 ms +/- 5    | 1.05x (4.6%) | 44 ms          |
| mpox/clade-ii/20  | 4.36 s +/- 0.04 | 4.23 s +/- 0.06 | 1.03x (2.9%) | 142 ms         |
| mpox/clade-ii/100 | 23.0 s +/- 0.1  | 22.4 s +/- 0.2  | 1.03x (2.9%) | 712 ms         |

The improvement is consistent but modest. `edge_subs` accounts for only ~3% of total ancestral reconstruction time; the optimization reduces its cost but cannot affect the dominant 36% math primitives or 13% GTR inference.

The absolute CPU saving scales with `edges * alignment_length`: mpox/clade-ii/100 saves 712ms from 38.8M fewer `ArrayView` row extractions and per-element comparisons.

## 5. Other optimization opportunities

### 5.1 Bitset for `range_contains`

Not profiled (perf_event_paranoid=4 on this host prevented CPU profiling). This is a structural observation, not a measured bottleneck.

`range_contains` (`range.rs:9`) performs a linear scan of `non_char` range tuples at each call. After commit `6396f86a`, it only runs at positions where parent and child disagree AND both states are canonical (<5% of positions on typical edges), so the total call count is low. For the common case of viral datasets with few contiguous gap/unknown ranges (R < 10), the linear scan is fast regardless.

A bitset would replace O(R) per lookup with O(1), but without profiling data showing `range_contains` as a measurable contributor, the benefit is speculative. Would need profiling on a dataset with many scattered gaps (TB, bacterial genomes with fragmented alignments) to justify.

### 5.2 `prof2seq` reuse (remove redundant argmax passes)

`assign_sequence()` calls `prof2seq()` which iterates the same profile matrix with `argmax_first` per row - the same work `edge_subs` does. When both functions run on the same node in the same pipeline (ancestral reconstruction output), the argmax computation is duplicated.

From the [OpenBLAS profiling report](dense-openblas-profiling.md), `argmax_first` accounts for 12.6% of total CPU across all call sites (marginal passes, `prof2seq`, `assign_sequence`, GTR inference). The `edge_subs` and `assign_sequence` call sites are a subset of that 12.6%. Caching the argmax index vector (or the derived sequence) on `DenseNodePartition` after the forward pass would eliminate the redundant passes.

For mpox/clade-ii/100 (197 nodes, 197k positions each), this removes ~39M redundant comparisons. Needs targeted profiling to measure the actual fraction of the 12.6% attributable to the duplicated calls.

### 5.3 Fused `edge_subs` + `edge_effective_length`

Structural observation, not profiled. `initial_guess_mixed()` (`optimize_unified.rs:748`) calls both `edge_subs()` and `edge_effective_length()` per edge. Both look up the same `non_char` ranges; `edge_effective_length` also clones and unions them via `range_union`. A fused function could compute both sub count and effective length in a single pass. Only applies to the `optimize` command path, not ancestral-only.

### 5.4 Column compression for dense mode

The highest-impact optimization for dense mode overall: group identical alignment columns and weight by multiplicity, matching what v0 does. For viral datasets, unique columns are typically 5-20% of total alignment length. This would reduce marginal pass cost, GTR inference, argmax, and `edge_subs` proportionally.

This is documented in existing known issues ([N-optimize-dense-iteration-slow](../port-known-issues/N-optimize-dense-iteration-slow.md)) and the [OpenBLAS profiling report](dense-openblas-profiling.md). It is the primary lever for dense mode performance improvement and would subsume the per-function micro-optimizations in 5.1-5.3.

## 6. Assessment

The vectorized argmax experiment (section 3) shows a consistent 3-5% wall-time improvement. The gain is inherently bounded by `edge_subs` being ~3% of total ancestral reconstruction time. Not worth the added complexity.

Column compression (5.4) is the structural fix that would reduce all O(L) costs proportionally. Argmax caching (5.2) is the most impactful incremental change within the current architecture.
