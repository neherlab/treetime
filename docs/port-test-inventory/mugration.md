# Mugration Test Coverage

[Back to index](_index.md)

## Golden Master (v0 parity)

| Test                                      | Datasets                                   | Status                   |
| :---------------------------------------- | :----------------------------------------- | :----------------------- |
| `test_gm_mugration_outputs`               | zika, zika_weights                         | passing                  |
| `test_gm_mugration_outputs_v1_divergence` | lassa, dengue, tb, rsv, mpox               | ignored (D1/D2)          |
| `test_gm_mugration_confidence_zika`       | zika                                       | passing (2e-2 tolerance) |
| `test_gm_mugration_confidence_outputs`    | zika_weights, lassa, dengue, tb, rsv, mpox | ignored (D1/D2)          |

Files: [packages/treetime/src/commands/mugration/**tests**/test_gm_mugration.rs](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs)

## Structural / Unit

| Test                                              | What it validates                                             |
| :------------------------------------------------ | :------------------------------------------------------------ |
| `test_execute_mugration_simple_tree`              | Basic 2-state tree: states, assignments, confidence structure |
| `test_execute_mugration_with_weights`             | Weight-based pi computation, weighted reconstruction          |
| `test_execute_mugration_with_missing_data`        | Missing trait handling                                        |
| `test_execute_mugration_with_pseudo_counts`       | Pseudo-count effect on pi with weights                        |
| `test_execute_mugration_sampling_bias_correction` | mu scaled by correction factor                                |
| `test_execute_mugration_rejects_single_state`     | Error on < 2 states                                           |

Files: [packages/treetime/src/commands/mugration/**tests**/test_run.rs](../../packages/treetime/src/commands/mugration/__tests__/test_run.rs)

## Algorithm Invariants

| Test                                           | What it validates                                  |
| :--------------------------------------------- | :------------------------------------------------- |
| `test_iterative_refinement_changes_model`      | mu/pi differ between iterations=0 and iterations=5 |
| `test_iterative_refinement_pi_reflects_data`   | pi shifts toward observed frequencies with free pi |
| `test_zero_iterations_preserves_initial_model` | Basic parameter validity at iterations=0           |

Files: [packages/treetime/src/commands/mugration/**tests**/test_run.rs](../../packages/treetime/src/commands/mugration/__tests__/test_run.rs)

## Brent Optimizer

| Test                                            | What it validates                  |
| :---------------------------------------------- | :--------------------------------- |
| `test_brent_minimize_quadratic`                 | Minimum of (x-3)^2                 |
| `test_brent_minimize_shifted_quadratic`         | Minimum of (x-0.7)^2+1             |
| `test_brent_minimize_cosine`                    | Minimum of cos(x) on [2,5]         |
| `test_brent_minimize_monotone_returns_boundary` | Monotone function returns boundary |

Files: [packages/treetime/src/commands/mugration/gtr_refinement.rs](../../packages/treetime/src/commands/mugration/gtr_refinement.rs)

## Partition / Discrete

| Test                           | What it validates          |
| :----------------------------- | :------------------------- |
| `test_new_partition`           | Constructor defaults       |
| `test_get_reconstructed_trait` | Argmax state extraction    |
| `test_get_confidence`          | Confidence profile access  |
| `test_get_log_lh`              | Log-likelihood access      |
| `test_argmax_first_1d`         | Deterministic tie-breaking |
| `test_normalize_*`             | Normalization edge cases   |

Files: [packages/treetime/src/representation/partition/discrete.rs](../../packages/treetime/src/representation/partition/discrete.rs)

## Other

| Test                       | What it validates                |
| :------------------------- | :------------------------------- |
| `test_comment_output_*`    | Nexus annotation format          |
| `test_discrete_marginal_*` | Forward-backward message passing |

Files: [packages/treetime/src/commands/mugration/**tests**/test_comment_output.rs](../../packages/treetime/src/commands/mugration/__tests__/test_comment_output.rs), [packages/treetime/src/commands/mugration/**tests**/test_discrete_marginal.rs](../../packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs)
