# Sparse substitution accessors

Production call sites for `fitch_subs` and `marginal_subs` on `SparseEdgePartition`, excluding definitions and tests.

## Definitions

All accessors defined on [`SparseEdgePartition`](../../packages/treetime/src/representation/payload/sparse.rs).

### Fitch subs

- [`fitch_subs(&self) -> &[Sub]`](../../packages/treetime/src/representation/payload/sparse.rs#L120)
- [`set_fitch_subs(&mut self, subs)`](../../packages/treetime/src/representation/payload/sparse.rs#L124)
- [`extend_fitch_subs(&mut self, subs)`](../../packages/treetime/src/representation/payload/sparse.rs#L129)
- [`invert_fitch_subs(&mut self)`](../../packages/treetime/src/representation/payload/sparse.rs#L134)
- [`chain_fitch_subs(&self, suffix) -> Result<Vec<Sub>>`](../../packages/treetime/src/representation/payload/sparse.rs#L141)
- [`with_fitch_subs(subs) -> Self`](../../packages/treetime/src/representation/payload/sparse.rs#L105) (constructor)
- [`with_fitch_subs_and_indels(subs, indels) -> Self`](../../packages/treetime/src/representation/payload/sparse.rs#L112) (constructor)

### Marginal subs

- [`marginal_subs(&self) -> Option<&[Sub]>`](../../packages/treetime/src/representation/payload/sparse.rs#L145)
- [`set_marginal_subs(&mut self, subs)`](../../packages/treetime/src/representation/payload/sparse.rs#L149)
- [`clear_marginal_subs(&mut self)`](../../packages/treetime/src/representation/payload/sparse.rs#L153)

## Fitch subs call sites (22)

### Ancestral reconstruction

| Accessor                          | Location                                                                      |
| :-------------------------------- | :---------------------------------------------------------------------------- |
| `extend_fitch_subs(subs)`         | [fitch.rs#L453](../../packages/treetime/src/commands/ancestral/fitch.rs#L453) |
| `fitch_subs()` iterate for output | [fitch.rs#L580](../../packages/treetime/src/commands/ancestral/fitch.rs#L580) |

### Marginal passes

| Accessor                                             | Location                                                                                                |
| :--------------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| `fitch_subs()` in `compute_marginal_subs_for_edge()` | [marginal_passes.rs#L100](../../packages/treetime/src/representation/partition/marginal_passes.rs#L100) |
| `fitch_subs()` forward pass                          | [marginal_passes.rs#L174](../../packages/treetime/src/representation/partition/marginal_passes.rs#L174) |
| `fitch_subs()` backward pass                         | [marginal_passes.rs#L272](../../packages/treetime/src/representation/partition/marginal_passes.rs#L272) |
| `fitch_subs()` child edge                            | [marginal_passes.rs#L364](../../packages/treetime/src/representation/partition/marginal_passes.rs#L364) |

### Sparse reroot

| Accessor                                     | Location                                                                                                |
| :------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| `fitch_subs()` read                          | [marginal_sparse.rs#L98](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L98)   |
| `fitch_subs()` parent edge                   | [marginal_sparse.rs#L171](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L171) |
| `chain_fitch_subs(child.fitch_subs())` merge | [marginal_sparse.rs#L189](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L189) |
| `set_fitch_subs(merged)` store               | [marginal_sparse.rs#L195](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L195) |
| `invert_fitch_subs()` reroot                 | [marginal_sparse.rs#L209](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L209) |
| `fitch_subs()` post-reroot                   | [marginal_sparse.rs#L235](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L235) |

### Prune

| Accessor                            | Location                                                              |
| :---------------------------------- | :-------------------------------------------------------------------- |
| `fitch_subs().len()` count          | [run.rs#L167](../../packages/treetime/src/commands/prune/run.rs#L167) |
| `fitch_subs()` pair A               | [run.rs#L416](../../packages/treetime/src/commands/prune/run.rs#L416) |
| `fitch_subs()` pair B               | [run.rs#L417](../../packages/treetime/src/commands/prune/run.rs#L417) |
| `fitch_subs()` pair A (alt path)    | [run.rs#L475](../../packages/treetime/src/commands/prune/run.rs#L475) |
| `fitch_subs()` pair B (alt path)    | [run.rs#L476](../../packages/treetime/src/commands/prune/run.rs#L476) |
| `set_fitch_subs()` shared on parent | [run.rs#L569](../../packages/treetime/src/commands/prune/run.rs#L569) |
| `set_fitch_subs()` remaining A      | [run.rs#L573](../../packages/treetime/src/commands/prune/run.rs#L573) |
| `set_fitch_subs()` remaining B      | [run.rs#L577](../../packages/treetime/src/commands/prune/run.rs#L577) |

### Topology cleanup

| Accessor                                                 | Location                                                                                            |
| :------------------------------------------------------- | :-------------------------------------------------------------------------------------------------- |
| `chain_fitch_subs(child.fitch_subs())` merge on collapse | [collapse.rs#L60](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L60) |
| `set_fitch_subs(merged)` store                           | [collapse.rs#L61](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L61) |

### Optimize

| Accessor                                      | Location                                                                                       |
| :-------------------------------------------- | :--------------------------------------------------------------------------------------------- |
| `fitch_subs().is_empty()` check               | [run.rs#L583](../../packages/treetime/src/commands/optimize/run.rs#L583)                       |
| `fitch_subs().iter().map(Sub::pos)` positions | [optimize_sparse.rs#L58](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L58) |
| `fitch_subs().iter().find()` lookup           | [optimize_sparse.rs#L66](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L66) |

### GTR inference

| Accessor                            | Location                                                                 |
| :---------------------------------- | :----------------------------------------------------------------------- |
| `fitch_subs()` parameter estimation | [sparse.rs#L75](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L75) |

### Output serialization

| Accessor                        | Location                                                                                  |
| :------------------------------ | :---------------------------------------------------------------------------------------- |
| `with_fitch_subs(subs.clone())` | [timetree.rs#L348](../../packages/treetime/src/representation/payload/timetree.rs#L348)   |
| `with_fitch_subs(subs.clone())` | [ancestral.rs#L246](../../packages/treetime/src/representation/payload/ancestral.rs#L246) |

## Marginal subs call sites (6)

### Marginal passes

| Accessor                           | Location                                                                                                |
| :--------------------------------- | :------------------------------------------------------------------------------------------------------ |
| `compute_marginal_subs_for_edge()` | [marginal_passes.rs#L81](../../packages/treetime/src/representation/partition/marginal_passes.rs#L81)   |
| `set_marginal_subs(marginal_subs)` | [marginal_passes.rs#L341](../../packages/treetime/src/representation/partition/marginal_passes.rs#L341) |

### Sparse reroot

| Accessor                     | Location                                                                                                |
| :--------------------------- | :------------------------------------------------------------------------------------------------------ |
| `clear_marginal_subs()`      | [marginal_sparse.rs#L146](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L146) |
| `clear_marginal_subs()`      | [marginal_sparse.rs#L212](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L212) |
| `marginal_subs()` match read | [marginal_sparse.rs#L269](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L269) |

### Output serialization

| Accessor                          | Location                                                                                  |
| :-------------------------------- | :---------------------------------------------------------------------------------------- |
| `set_marginal_subs(subs.clone())` | [timetree.rs#L349](../../packages/treetime/src/representation/payload/timetree.rs#L349)   |
| `set_marginal_subs(subs.clone())` | [ancestral.rs#L247](../../packages/treetime/src/representation/payload/ancestral.rs#L247) |
