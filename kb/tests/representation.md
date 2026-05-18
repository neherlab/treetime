# Representation Tests

[Back to index](README.md)

## Summary

| Category                                                             | Type                 |
| -------------------------------------------------------------------- | -------------------- |
| [Substitution composition](#substitution-composition)                | Unit + Parameterized |
| [Discrete states](#discrete-states)                                  | Unit                 |
| [Normalize from log (dense 2D)](#normalize-from-log-dense-2d)        | Unit                 |
| [Normalize inplace (dense 2D)](#normalize-inplace-dense-2d)          | Unit                 |
| [Normalize (discrete 1D)](#normalize-discrete-1d)                    | Unit                 |
| [Per-site rate propagation](#per-site-rate-propagation)              | Unit                 |
| [Topology cleanup / edge collapse](#topology-cleanup--edge-collapse) | Unit                 |
| [Topology cleanup / reroot](#topology-cleanup--reroot)               | Unit                 |
| [Payload: ancestral annotation](#payload-ancestral-annotation)       | Unit                 |
| [Payload: discrete data](#payload-discrete-data)                     | Unit                 |
| [Payload: timetree annotation](#payload-timetree-annotation)         | Unit                 |

Support files (helpers only, no tests): [`packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_collapse_edge.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_collapse_edge.rs) contains inline helpers `c()`, `sub()`, `populate_test_nodes()`, `make_sparse_partition()`.

---

## Substitution Composition

**Test:** [`packages/treetime/src/partition/__tests__/test_partition_marginal_sparse.rs`](../../packages/treetime/src/partition/__tests__/test_partition_marginal_sparse.rs)

**Impl:** [`packages/treetime/src/seq/mutation.rs`](../../packages/treetime/src/seq/mutation.rs)

| Test                                                   | Purpose                                     |
| ------------------------------------------------------ | ------------------------------------------- |
| `test_compose_substitutions_empty_both`                | Empty parent and child yields empty         |
| `test_compose_substitutions_empty_parent`              | Empty parent preserves child subs           |
| `test_compose_substitutions_empty_child`               | Empty child preserves parent subs           |
| `test_compose_substitutions_no_overlap`                | Non-overlapping positions merged            |
| `test_compose_substitutions_chain`                     | Same-position subs compose (A->G->T = A->T) |
| `test_compose_substitutions_cancellation`              | Reverse sub cancels (A->G->A = none)        |
| `test_compose_substitutions_mixed`                     | Chain, keep, add, cancel combined           |
| `test_compose_substitutions_single_position` (4 cases) | Single-position chain and cancel variants   |

**Algorithm:** Two-pointer merge composition of parent and child substitution lists by position: chain, cancellation, or passthrough for non-overlapping positions.

---

## Discrete States

**Test:** [`packages/treetime/src/partition/discrete_states.rs`](../../packages/treetime/src/partition/discrete_states.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/discrete_states.rs`](../../packages/treetime/src/partition/discrete_states.rs)

| Test                                            | Purpose                                                  |
| ----------------------------------------------- | -------------------------------------------------------- |
| `test_from_values_sorts_and_deduplicates`       | from_values sorts, deduplicates, excludes missing marker |
| `test_get_index_returns_none_for_missing`       | Missing marker returns None index                        |
| `test_get_index_returns_correct_index`          | Indices match sorted order                               |
| `test_get_name_returns_correct_name`            | Names match sorted order                                 |
| `test_is_missing`                               | Missing marker detected, non-missing values not flagged  |
| `test_get_index_returns_none_for_unknown_value` | Unknown value returns None                               |
| `test_empty_values`                             | Empty input produces empty DiscreteStates                |
| `test_missing_marker`                           | Custom missing marker stored and returned                |

---

## Normalize from Log (Dense 2D)

**Test:** [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs)

| Test                                                      | Purpose                                                            |
| --------------------------------------------------------- | ------------------------------------------------------------------ |
| `test_normalize_from_log_equal_probs`                     | Equal log-probs produce uniform distribution, finite log_lh        |
| `test_normalize_from_log_descending`                      | Unequal log-probs: analytically verified output and norm           |
| `test_normalize_from_log_all_neg_inf_single_row`          | All-NEG_INFINITY single row: uniform 1/4, log_lh=-inf              |
| `test_normalize_from_log_all_neg_inf_multiple_rows`       | All-NEG_INFINITY across 3 rows: uniform fallback each              |
| `test_normalize_from_log_mixed_neg_inf_and_finite_rows`   | Degenerate row 0 + normal row 1: row-level fallback, log_lh=-inf   |
| `test_normalize_from_log_mixed_finite_neg_inf_within_row` | Mixed finite/-inf within a row: only finite states get probability |
| `test_normalize_from_log_large_negative_values`           | Numerical stability: offset -1000 matches reference                |
| `test_normalize_from_log_three_states`                    | All-NEG_INFINITY fallback with 3 states: uniform 1/3               |

**Algorithm:** 2D logsumexp normalization matching v0's `normalize_profile(log=True)`. Guards all-`-inf` rows with uniform fallback, consistent with sparse path's `logsumexp_normalize`.

---

## Normalize Inplace (Dense 2D)

**Test:** [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs) (inline `#[cfg(test)]`, same module)

**Impl:** [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs)

| Test                                                | Purpose                                              |
| --------------------------------------------------- | ---------------------------------------------------- |
| `test_normalize_inplace_normal_rows`                | Normal rows: correct normalization and LH            |
| `test_normalize_inplace_zero_row_returns_uniform`   | All-zero row: uniform fallback, log_lh=-inf          |
| `test_normalize_inplace_mixed_zero_and_normal_rows` | Zero row + normal row: per-row handling, log_lh=-inf |
| `test_normalize_inplace_nan_row_returns_uniform`    | NaN row: uniform fallback, log_lh=-inf               |
| `test_normalize_inplace_inf_row_returns_uniform`    | Inf row: uniform fallback, log_lh=-inf               |

**Algorithm:** `normalize_inplace` guards zero-sum and non-finite rows with the same uniform fallback as `normalize_from_log`.

---

## Normalize (Discrete 1D)

**Test:** [`packages/treetime/src/partition/discrete.rs`](../../packages/treetime/src/partition/discrete.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/discrete.rs`](../../packages/treetime/src/partition/discrete.rs)

| Test                                                   | Purpose                                     |
| ------------------------------------------------------ | ------------------------------------------- |
| `test_new_partition`                                   | Constructor defaults                        |
| `test_get_reconstructed_trait`                         | Argmax state extraction                     |
| `test_get_confidence`                                  | Confidence profile access                   |
| `test_get_log_lh`                                      | Log-likelihood access                       |
| `test_discrete_argmax_first`                           | Deterministic tie-breaking (shared utility) |
| `test_normalize_inplace_1d_zero_sum_returns_error`     | Zero-sum normalization returns error        |
| `test_normalize_from_log_1d_all_neg_inf_returns_error` | All-neg-inf log normalization returns error |
| `test_normalize_inplace_1d_valid_input`                | Valid normalization sums to 1               |
| `test_normalize_from_log_1d_valid_input`               | Valid log normalization sums to 1           |

Also cross-referenced from [Mugration Tests](mugration.md#partition--discrete).

---

## Per-Site Rate Propagation

**Test:** [`packages/treetime/src/partition/marginal_helpers.rs`](../../packages/treetime/src/partition/marginal_helpers.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/marginal_helpers.rs`](../../packages/treetime/src/partition/marginal_helpers.rs)

| Test                                   | Purpose                                                                    |
| -------------------------------------- | -------------------------------------------------------------------------- |
| `test_propagate_raw_per_site_forward`  | Forward propagation with per-site rates matches individual expQt_with_rate |
| `test_propagate_raw_per_site_backward` | Backward propagation (transpose) with per-site rates matches individual    |

**Algorithm:** Per-site rate-scaled transition matrix propagation for Felsenstein pruning with site-rate variation.

---

## Topology Cleanup / Edge Collapse

**Test:** [`packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_collapse_edge.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_collapse_edge.rs)

**Impl:** [`packages/treetime/src/partition/algo/topology_cleanup/collapse.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/collapse.rs)

| Test                                                      | Purpose                                                             |
| --------------------------------------------------------- | ------------------------------------------------------------------- |
| `test_topology_collapse_edge_sparse_composes_subs`        | Substitutions composed correctly on sparse collapse                 |
| `test_topology_collapse_edge_dense_cleanup`               | Stale dense partition data removed after collapse                   |
| `test_topology_collapse_edge_branch_length_sum`           | Branch lengths summed correctly (non-zero collapsed + child)        |
| `test_topology_collapse_edge_branch_length_sum_with_zero` | Branch lengths summed correctly (zero collapsed + child unchanged)  |
| `test_topology_collapse_edge_indel_concatenation`         | Collapsed-edge indels prepended to child indels                     |
| `test_topology_collapse_edge_reversion_cancels`           | Forward + reverse substitution at same position cancel to no change |
| `test_topology_collapse_edge_no_partitions`               | Graph-only collapse with no partitions still rewires topology       |
| `test_topology_collapse_edge_multiple_sparse_partitions`  | Multiple sparse partitions each receive composed substitutions      |

Also cross-referenced from [Branch Optimization Tests](optimization.md#topology-cleanup-in-loop).

---

## Topology Cleanup / Reroot

**Test:** [`packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_reroot.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/__tests__/test_reroot.rs)

**Impl:** [`packages/treetime/src/partition/algo/topology_cleanup/reroot.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/reroot.rs)

| Test                                                          | Purpose                                                               |
| ------------------------------------------------------------- | --------------------------------------------------------------------- |
| `test_reroot_split_edge_divides_branch_length`                | Edge split distributes branch length at specified fractional position |
| `test_reroot_split_edge_at_midpoint`                          | Midpoint split produces equal-length child and parent edges           |
| `test_reroot_apply_reroot_topology_inverts_path`              | Single-hop reroot inverts one edge and moves root                     |
| `test_reroot_apply_reroot_topology_multi_hop`                 | Multi-hop reroot inverts two edges through internal node              |
| `test_reroot_apply_reroot_topology_preserves_leaf_count`      | Leaf count unchanged after topology reroot                            |
| `test_reroot_remove_node_if_trivial_merges_edges`             | Degree-2 node removed with summed branch lengths                      |
| `test_reroot_remove_node_if_trivial_non_trivial_returns_none` | Non-trivial nodes (root, bifurcating internal) return None            |
| `test_reroot_full_reroot_and_cleanup_preserves_topology`      | Full reroot + trivial removal preserves all leaves and taxa           |
| `test_reroot_remove_trivial_with_partial_branch_lengths`      | Partial branch lengths (Some + None) produce Some after merge         |
| `test_reroot_full_cycle_branch_length_conservation`           | Total branch length conserved across reroot + trivial node removal    |

---

## Payload: Ancestral Annotation

**Test:** [`packages/treetime/src/partition/payload/ancestral.rs`](../../packages/treetime/src/partition/payload/ancestral.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/payload/ancestral.rs`](../../packages/treetime/src/partition/payload/ancestral.rs)

| Test                                                       | Purpose                                              |
| ---------------------------------------------------------- | ---------------------------------------------------- |
| `test_annotate_branch_mutations_formats_1_based_positions` | Mutations formatted with 1-based positions           |
| `test_annotate_branch_mutations_empty_partitions`          | Empty partition list produces no mutations           |
| `test_annotate_branch_mutations_no_mutations_on_edge`      | Edge with no subs produces None                      |
| `test_annotate_branch_mutations_sorts_by_position`         | Multiple mutations sorted by position                |
| `test_annotate_branch_mutations_multi_partition_merge`     | Mutations from multiple partitions merged and sorted |

---

## Payload: Discrete Data

**Test:** [`packages/treetime/src/partition/payload/discrete.rs`](../../packages/treetime/src/partition/payload/discrete.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/payload/discrete.rs`](../../packages/treetime/src/partition/payload/discrete.rs)

| Test                                         | Purpose                                             |
| -------------------------------------------- | --------------------------------------------------- |
| `test_from_observed_creates_one_hot_profile` | Observed index creates one-hot profile              |
| `test_missing_creates_uniform_profile`       | Missing creates uniform profile                     |
| `test_default_node_data`                     | Default node data has empty profile and zero log_lh |
| `test_default_edge_data`                     | Default edge data has empty arrays and zero log_lh  |

---

## Payload: Timetree Annotation

**Test:** [`packages/treetime/src/partition/payload/timetree.rs`](../../packages/treetime/src/partition/payload/timetree.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/partition/payload/timetree.rs`](../../packages/treetime/src/partition/payload/timetree.rs)

| Test                                                           | Purpose                                                         |
| -------------------------------------------------------------- | --------------------------------------------------------------- |
| `test_timetree_annotate_branch_mutations_populates_base_field` | annotate_branch_mutations writes to NodeTimetree.base.mutations |
| `test_timetree_nwk_comments_include_mutations_and_date`        | nwk_comments includes both mutations and date annotations       |
| `test_timetree_nexus_output_includes_mutations_and_date`       | Nexus serialization carries mutations and date annotations      |

---

## Files Without Tests

Production source files in `representation/` with no inline or dedicated tests:

| File                           | Content                          |
| ------------------------------ | -------------------------------- |
| `partition/fitch.rs`           | Fitch partition state management |
| `partition/fitch_config.rs`    | Fitch configuration              |
| `partition/likelihood.rs`      | Likelihood computation helpers   |
| `partition/marginal_passes.rs` | Marginal forward/backward passes |
| `partition/marginal_sparse.rs` | Sparse marginal partition        |
| `partition/timetree.rs`        | Timetree partition               |
| `partition/traits.rs`          | Partition trait definitions      |
| `payload/dense.rs`             | Dense payload structs            |
| `payload/sparse.rs`            | Sparse payload structs           |
| `algo/infer_dense.rs`          | Dense inference algorithm        |

These are exercised indirectly through command-level tests (ancestral, optimize, timetree) but have no direct unit tests.

## Known Test Deficiencies

Specific functions and code paths identified as needing direct test coverage:

| Entity                              | File                                                          | Deficiency                                                                                                                      |
| ----------------------------------- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `combine_messages()`                | `partition/marginal_helpers.rs`                               | Core sparse marginal algorithm with no direct tests. Only indirect coverage via command integration tests.                      |
| `reconcile_topology()`              | `partition/marginal_sparse.rs`, `partition/marginal_dense.rs` | Adds/removes partition entries after graph topology edits with no tests asserting correctness.                                  |
| `propagate_raw_per_site()` tests    | `partition/marginal_helpers.rs`                               | Tests use circular oracle (same `expQt_with_rate()` as SUT). Need hand-computed expected values.                                |
| `collapse_edge()` composition tests | `algo/topology_cleanup/__tests__/test_collapse_edge.rs`       | Assert substitution count, not exact content. Composition semantics verified separately in `test_partition_marginal_sparse.rs`. |
