# Distribution Operations Tests

[Back to index](README.md)

## Summary

| Category                                                                | Crate                 | Type |
| ----------------------------------------------------------------------- | --------------------- | ---- |
| [Distribution core](#distribution-core)                                 | treetime-distribution | Unit |
| [Distribution ops - multiply](#distribution-operations---multiply)      | treetime-distribution | Unit |
| [Distribution ops - divide](#distribution-operations---divide)          | treetime-distribution | Unit |
| [Distribution ops - convolve](#distribution-operations---convolve)      | treetime-distribution | Unit |
| [Distribution ops - negation](#distribution-operations---negation)      | treetime-distribution | Unit |
| [Distribution ops - scalar](#distribution-operations---scalar-multiply) | treetime-distribution | Unit |
| [Distribution ops - bounds](#distribution-operations---time-bounds)     | treetime-distribution | Unit |
| [Scaled distributions](#scaled-distributions)                           | treetime-distribution | Unit |
| [Gaussian product](#root-level-tests)                                   | treetime-distribution | Unit |
| [Quantile](#root-level-tests)                                           | treetime-distribution | Unit |
| [Gaussian analytical](#treetime-analytical)                             | treetime-analytical   | Unit |
| [Exponential](#exponential-inline-tests)                                | treetime-analytical   | Unit |
| [Gaussian-exponential](#gaussian-exponential-inline-tests)              | treetime-analytical   | Unit |
| [Array multiplication](#treetime-ops)                                   | treetime-ops          | Unit |
| [Convolution](#convolution-inline-tests)                                | treetime-ops          | Unit |

---

## treetime-distribution

### Distribution Core

**Directory:** [`distribution_core/__tests__/`](../../packages/treetime-distribution/src/distribution_core/__tests__/)

**Test:** [`packages/treetime-distribution/src/distribution_core/__tests__/test_distribution.rs`](../../packages/treetime-distribution/src/distribution_core/__tests__/test_distribution.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_core/distribution.rs`](../../packages/treetime-distribution/src/distribution_core/distribution.rs)

| Test                                        | Purpose                        |
| ------------------------------------------- | ------------------------------ |
| `test_distribution_max_value`               | Max value across variants      |
| `test_distribution_scale_by`                | Scale amplitude by factor      |
| `test_distribution_normalize`               | Normalize to max=1             |
| `test_distribution_normalize_empty_on_zero` | Zero values normalize to Empty |

**Test:** [`packages/treetime-distribution/src/distribution_core/__tests__/test_function.rs`](../../packages/treetime-distribution/src/distribution_core/__tests__/test_function.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_core/function.rs`](../../packages/treetime-distribution/src/distribution_core/function.rs)

| Test                                          | Purpose                          |
| --------------------------------------------- | -------------------------------- |
| `test_distribution_function_interpolate`      | Linear interpolation at point    |
| `test_distribution_function_interpolate_many` | Interpolation at multiple points |

---

### Distribution Operations - Multiply

**Directory:** [`distribution_ops/__tests__/`](../../packages/treetime-distribution/src/distribution_ops/__tests__/)

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs)

| Test                                                            | Purpose                                      |
| --------------------------------------------------------------- | -------------------------------------------- |
| `test_multiply_formula_function_returns_function`               | Formula \* Function = Function               |
| `test_multiply_function_formula_commutative`                    | Function \* Formula commutativity            |
| `test_multiply_function_function_non_overlapping_returns_empty` | Non-overlapping Function \* Function = Empty |
| `test_multiply_function_function_overlapping_uses_intersection` | Overlapping uses intersection range          |
| `test_multiply_chain_with_normalize_no_underflow`               | Chain multiply + normalize resists underflow |
| `test_multiplication_eval_range`                                | Parameterized range intersection (7 cases)   |

---

### Distribution Operations - Divide

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/divide.rs`](../../packages/treetime-distribution/src/distribution_ops/divide.rs)

| Test                                               | Purpose                            |
| -------------------------------------------------- | ---------------------------------- |
| `test_divide_empty_by_any`                         | Empty / any = Empty                |
| `test_divide_by_empty_fails`                       | Division by empty errors           |
| `test_divide_point_by_function`                    | Point / Function at location       |
| `test_divide_function_by_function`                 | Function / Function elementwise    |
| `test_divide_by_zero_handled`                      | Zero divisor uses tiny number      |
| `test_divide_range_by_function_full_overlap`       | Range / Function full overlap      |
| `test_divide_range_by_function_partial_overlap`    | Range / Function partial overlap   |
| `test_divide_range_by_function_no_overlap`         | Range / Function no overlap        |
| `test_divide_function_by_function_same_grid`       | Same grid elementwise division     |
| `test_divide_function_by_function_different_grids` | Different grids with interpolation |

---

### Distribution Operations - Convolve

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs)

| Test                                                   | Purpose                          |
| ------------------------------------------------------ | -------------------------------- |
| `test_convolution_empty`                               | Empty convolution returns Empty  |
| `test_convolution_point_point`                         | Point + Point = shifted Point    |
| `test_convolution_range_range_triangle`                | Equal-width ranges make triangle |
| `test_convolution_range_range_trapezoid_non_uniform`   | Unequal ranges make trapezoid    |
| `test_convolution_range_range_trapezoid_uniform`       | Uniform trapezoid case           |
| `test_convolution_point_function`                      | Point shifts Function            |
| `test_convolution_range_function`                      | Range smooths Function           |
| `test_convolution_point_range`                         | Point + Range = shifted Range    |
| `test_convolution_range_point`                         | Range + Point = shifted Range    |
| `test_convolution_function_function_basic`             | Basic function-function convolve |
| `test_convolution_function_function_single_points`     | Single-point functions convolve  |
| `test_convolution_function_function_empty`             | Empty function returns Empty     |
| `test_convolution_function_function_different_spacing` | Different grid spacing           |
| `test_convolution_function_function_zero_width`        | Zero-width function convolve     |
| `test_backward_pass_temporal_direction`                | Backward pass time direction     |
| `test_forward_pass_temporal_direction`                 | Forward pass time direction      |
| `test_convolution_with_uncertainty`                    | Uncertainty distributions        |
| `test_convolution_convolve_small_dx_function_function` | Small dx regression test         |
| `test_convolution_convolve_small_dx_range_function`    | Small dx range regression test   |

---

### Distribution Operations - Negation

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_negation.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_negation.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/negate.rs`](../../packages/treetime-distribution/src/distribution_ops/negate.rs)

| Test                           | Purpose                          |
| ------------------------------ | -------------------------------- |
| `test_negate_empty`            | Negating Empty stays Empty       |
| `test_negate_point`            | Point location negated           |
| `test_negate_point_zero`       | Zero point unchanged             |
| `test_negate_range`            | Range bounds negated and swapped |
| `test_negate_range_symmetric`  | Symmetric range unchanged        |
| `test_negate_function`         | Function grid reversed           |
| `test_negate_inplace_point`    | In-place point negation          |
| `test_negate_inplace_range`    | In-place range negation          |
| `test_negate_inplace_function` | In-place function negation       |

---

### Distribution Operations - Scalar Multiply

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_scalar_multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_scalar_multiply.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/scalar_multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/scalar_multiply.rs)

| Test                                                                 | Purpose                        |
| -------------------------------------------------------------------- | ------------------------------ |
| `test_distribution_scalar_multiplication_function_positive_scalar`   | Function \* positive scalar    |
| `test_distribution_scalar_multiplication_function_zero_scalar`       | Function \* zero = zeros       |
| `test_distribution_scalar_multiplication_function_negative_scalar`   | Function \* negative scalar    |
| `test_distribution_scalar_multiplication_point_positive_scalar`      | Point \* positive scalar       |
| `test_distribution_scalar_multiplication_point_zero_scalar`          | Point \* zero = zero amplitude |
| `test_distribution_scalar_multiplication_point_negative_scalar`      | Point \* negative scalar       |
| `test_distribution_scalar_multiplication_range_positive_scalar`      | Range \* positive scalar       |
| `test_distribution_scalar_multiplication_range_zero_scalar`          | Range \* zero = zero amplitude |
| `test_distribution_scalar_multiplication_range_negative_scalar`      | Range \* negative scalar       |
| `test_distribution_scalar_multiplication_empty`                      | Empty \* scalar = Empty        |
| `test_distribution_scalar_multiplication_function_fractional_scalar` | Function \* 0.1 fractional     |

---

### Distribution Operations - Time Bounds

**Test:** [`packages/treetime-distribution/src/distribution_ops/__tests__/test_time_bounds.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_time_bounds.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_ops/time_bounds.rs`](../../packages/treetime-distribution/src/distribution_ops/time_bounds.rs)

**Parameterized: `test_distribution_time_bounds_union`**

| Case                          | Purpose                        |
| ----------------------------- | ------------------------------ |
| `two_points_ordered`          | Union of ordered points        |
| `two_points_reversed`         | Union of reversed points       |
| `overlapping_ranges`          | Overlapping ranges union       |
| `overlapping_ranges_reversed` | Reversed overlapping ranges    |
| `disjoint_ranges`             | Disjoint ranges span gap       |
| `nested_ranges`               | Nested ranges use outer bounds |
| `point_and_range`             | Point and range union          |
| `two_functions`               | Function domains union         |
| `same_point`                  | Same point identity            |
| `negative_and_positive`       | Negative and positive ranges   |

**Parameterized: `test_distribution_time_bounds_intersection`**

| Case                          | Purpose                         |
| ----------------------------- | ------------------------------- |
| `overlapping_ranges`          | Overlapping ranges intersection |
| `overlapping_ranges_reversed` | Reversed overlap same result    |
| `disjoint_ranges`             | Disjoint returns None           |
| `disjoint_ranges_reversed`    | Reversed disjoint returns None  |
| `nested_outer_first`          | Nested returns inner bounds     |
| `nested_inner_first`          | Reversed nesting same result    |
| `point_inside_range`          | Point in range = point bounds   |
| `range_contains_point`        | Range with point = point bounds |
| `point_outside_range`         | Point outside returns None      |
| `two_functions`               | Function domain intersection    |
| `adjacent_ranges`             | Adjacent ranges touch at point  |
| `same_point`                  | Same point = point bounds       |
| `different_points`            | Different points returns None   |

**Parameterized: `test_distribution_time_bounds_contains`**

| Case                         | Purpose                        |
| ---------------------------- | ------------------------------ |
| `outer_contains_inner`       | Outer range contains inner     |
| `inner_not_contains_outer`   | Inner does not contain outer   |
| `range_contains_point`       | Range contains point inside    |
| `point_not_contains_range`   | Point does not contain range   |
| `overlapping_not_contained`  | Overlap without containment    |
| `same_range`                 | Same range contains itself     |
| `same_point`                 | Same point contains itself     |
| `inner_extends_left`         | Left extension not contained   |
| `inner_extends_right`        | Right extension not contained  |
| `function_contains_function` | Wider function contains narrow |

**Parameterized: `test_distribution_time_bounds_overlaps`**

| Case                  | Purpose                          |
| --------------------- | -------------------------------- |
| `overlapping_ranges`  | Overlapping ranges overlap       |
| `disjoint_ranges`     | Disjoint ranges no overlap       |
| `adjacent_ranges`     | Adjacent ranges overlap at point |
| `nested_ranges`       | Nested ranges overlap            |
| `point_inside_range`  | Point inside range overlaps      |
| `point_outside_range` | Point outside no overlap         |
| `same_point`          | Same point overlaps              |
| `different_points`    | Different points no overlap      |

**Plain tests:**

| Test                                                           | Purpose                          |
| -------------------------------------------------------------- | -------------------------------- |
| `test_distribution_time_bounds_union_commutativity`            | Union is commutative             |
| `test_distribution_time_bounds_intersection_commutativity`     | Intersection is commutative      |
| `test_distribution_time_bounds_overlaps_symmetry`              | Overlaps is symmetric            |
| `test_distribution_time_bounds_union_associativity`            | Union is associative             |
| `test_distribution_time_bounds_intersection_associativity`     | Intersection is associative      |
| `test_distribution_time_bounds_idempotence`                    | Union/intersection idempotent    |
| `test_distribution_time_bounds_contains_reflexivity`           | Contains is reflexive            |
| `test_distribution_time_bounds_contains_transitivity`          | Contains is transitive           |
| `test_distribution_time_bounds_overlaps_reflexivity`           | Overlaps is reflexive            |
| `test_distribution_time_bounds_intersection_implies_overlap`   | Intersection implies overlap     |
| `test_distribution_time_bounds_contains_implies_overlap`       | Contains implies overlap         |
| `test_distribution_time_bounds_union_contains_both`            | Union contains both inputs       |
| `test_distribution_time_bounds_intersection_contained_by_both` | Intersection contained by both   |
| `test_distribution_time_bounds_negative_ranges`                | Negative range arithmetic        |
| `test_distribution_time_bounds_zero_width_intervals`           | Zero-width interval operations   |
| `test_distribution_time_bounds_single_point_overlap`           | Single-point overlap at boundary |

---

### Scaled Distributions

**Directory:** [`distribution_scaled/__tests__/`](../../packages/treetime-distribution/src/distribution_scaled/__tests__/)

**Test:** [`packages/treetime-distribution/src/distribution_scaled/__tests__/test_distribution_scaled.rs`](../../packages/treetime-distribution/src/distribution_scaled/__tests__/test_distribution_scaled.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_scaled/scaled.rs`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs)

| Test                                             | Purpose                             |
| ------------------------------------------------ | ----------------------------------- |
| `test_scaled_distribution_from_plain_round_trip` | Plain -> Scaled -> Plain round trip |
| `test_scaled_distribution_log_scale`             | Log scale tracks amplitude          |
| `test_scaled_distribution_inner_is_normalized`   | Inner distribution max = 1          |
| `test_scaled_distribution_empty`                 | Empty distribution handling         |
| `test_scaled_distribution_default`               | Default is empty                    |
| `test_scaled_distribution_renormalize`           | Renormalize preserves scale         |
| `test_scaled_distribution_from_parts`            | Construct from parts                |

**Test:** [`packages/treetime-distribution/src/distribution_scaled/__tests__/test_multiply.rs`](../../packages/treetime-distribution/src/distribution_scaled/__tests__/test_multiply.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_scaled/multiply.rs`](../../packages/treetime-distribution/src/distribution_scaled/multiply.rs)

| Test                                                                    | Purpose                          |
| ----------------------------------------------------------------------- | -------------------------------- |
| `test_scaled_distribution_multiply_points_same_location`                | Point \* Point same location     |
| `test_scaled_distribution_multiply_points_different_location`           | Point \* Point different = Empty |
| `test_scaled_distribution_multiply_functions`                           | Function \* Function elementwise |
| `test_scaled_distribution_multiply_many_no_underflow`                   | 100 small points no underflow    |
| `test_scaled_distribution_multiply_empty_propagates`                    | Empty \* any = Empty             |
| `test_scaled_distribution_multiply_preserves_normalization`             | Product inner max = 1            |
| `test_scaled_distribution_multiply_many_empty`                          | Multiply-many empty = Empty      |
| `test_scaled_distribution_multiply_many_single`                         | Multiply-many single = identity  |
| `test_scaled_distribution_multiply_many_functions_same_grid`            | Three functions same grid        |
| `test_scaled_distribution_multiply_many_functions_underflow_resistance` | 10 tiny functions no underflow   |
| `test_scaled_distribution_multiply_many_mixed_types_fallback`           | Mixed types use fallback path    |
| `test_scaled_distribution_multiply_many_different_grids_fallback`       | Different grids use fallback     |

**Test:** [`packages/treetime-distribution/src/distribution_scaled/__tests__/test_divide.rs`](../../packages/treetime-distribution/src/distribution_scaled/__tests__/test_divide.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_scaled/divide.rs`](../../packages/treetime-distribution/src/distribution_scaled/divide.rs)

| Test                                                      | Purpose                      |
| --------------------------------------------------------- | ---------------------------- |
| `test_scaled_distribution_divide_empty_by_any`            | Empty / any = Empty          |
| `test_scaled_distribution_divide_by_empty_fails`          | Division by empty errors     |
| `test_scaled_distribution_divide_functions`               | Function / Function quotient |
| `test_scaled_distribution_divide_inverse_of_multiply`     | (a\*b)/b recovers a          |
| `test_scaled_distribution_divide_preserves_normalization` | Quotient inner max = 1       |

**Test:** [`packages/treetime-distribution/src/distribution_scaled/__tests__/test_convolve.rs`](../../packages/treetime-distribution/src/distribution_scaled/__tests__/test_convolve.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_scaled/convolve.rs`](../../packages/treetime-distribution/src/distribution_scaled/convolve.rs)

| Test                                                        | Purpose                         |
| ----------------------------------------------------------- | ------------------------------- |
| `test_scaled_distribution_convolve_points`                  | Point + Point convolution       |
| `test_scaled_distribution_convolve_functions`               | Function + Function convolution |
| `test_scaled_distribution_convolve_empty`                   | Empty convolution = Empty       |
| `test_scaled_distribution_convolve_preserves_normalization` | Convolution inner max = 1       |
| `test_scaled_distribution_convolve_symmetric`               | Convolution commutativity       |

---

### Root-Level Tests

**Test:** [`packages/treetime-distribution/src/__tests__/test_gaussian_product.rs`](../../packages/treetime-distribution/src/__tests__/test_gaussian_product.rs)

**Impl:**

- [`packages/treetime-distribution/src/distribution_scaled/scaled.rs`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs)
- [`packages/treetime-distribution/src/distribution_scaled/multiply.rs`](../../packages/treetime-distribution/src/distribution_scaled/multiply.rs)
- [`packages/treetime-distribution/src/distribution_scaled/convolve.rs`](../../packages/treetime-distribution/src/distribution_scaled/convolve.rs)

| Test                                               | Purpose                          |
| -------------------------------------------------- | -------------------------------- |
| `test_gaussian_scaled_distribution_round_trip`     | Plain -> Scaled round trip       |
| `test_gaussian_point_multiplication_same_location` | Point \* Point same location     |
| `test_gaussian_underflow_resistance_points`        | 100 small points no underflow    |
| `test_gaussian_product_two_identical`              | Two N(0,1) product width         |
| `test_gaussian_product_ten`                        | Ten N(0,1) product width         |
| `test_gaussian_underflow_resistance_functions`     | 50 small functions no underflow  |
| `test_gaussian_convolution_basic`                  | Two N(0,1) convolution width     |
| `test_gaussian_analytical_formula`                 | Analytical product formula       |
| `test_gaussian_product_matches_analytical`         | Numerical matches analytical     |
| `test_gaussian_normalization_preserved`            | Amplitude tracking via log_scale |

**Test:** [`packages/treetime-distribution/src/__tests__/test_quantile.rs`](../../packages/treetime-distribution/src/__tests__/test_quantile.rs)

**Impl:** [`packages/treetime-distribution/src/distribution_core/distribution.rs`](../../packages/treetime-distribution/src/distribution_core/distribution.rs)

| Test                                          | Purpose                         |
| --------------------------------------------- | ------------------------------- |
| `test_quantile_invalid_p_returns_none`        | Invalid p returns None          |
| `test_quantile_empty_returns_none`            | Empty distribution returns None |
| `test_quantile_point_returns_point_location`  | Point quantile = point location |
| `test_quantile_range_linear_interpolation`    | Range quantile interpolation    |
| `test_quantile_uniform_function`              | Uniform function quantiles      |
| `test_quantile_single_point_function`         | Single-point function quantile  |
| `test_quantile_function_median`               | Triangular distribution median  |
| `test_quantile_function_boundaries`           | p=0 and p=1 boundary values     |
| `test_confidence_interval_95_percent`         | 95% confidence interval         |
| `test_confidence_interval_empty_returns_none` | Empty CI returns None           |

---

## treetime-analytical

**Directory:** [`__tests__/`](../../packages/treetime-analytical/src/__tests__/)

**Test:** [`packages/treetime-analytical/src/__tests__/gaussian.rs`](../../packages/treetime-analytical/src/__tests__/gaussian.rs)

**Impl:** [`packages/treetime-analytical/src/gaussian.rs`](../../packages/treetime-analytical/src/gaussian.rs)

| Test                                         | Purpose                        |
| -------------------------------------------- | ------------------------------ |
| `test_gaussian_product_params_single`        | Single Gaussian product params |
| `test_gaussian_product_params_two_identical` | Two identical product params   |
| `test_gaussian_product_params_two_shifted`   | Two shifted product params     |
| `test_gaussian_product_params_empty`         | Empty product params           |
| `test_gaussian_product_shape`                | Product shape on grid          |
| `test_gaussian_evaluate`                     | Gaussian evaluation at points  |
| `test_gaussian_convolution_same_width`       | Same-width convolution peak    |
| `test_gaussian_pdf_at_mean`                  | PDF value at mean              |
| `test_gaussian_pdf_symmetry`                 | PDF symmetry around mean       |
| `test_gaussian_pdf_grid`                     | PDF evaluation on grid         |
| `test_gaussian_convolution_pdf_centered`     | Convolution PDF at center      |
| `test_gaussian_convolution_pdf_shifted`      | Convolution PDF shifted        |

**Note:** `test_gaussian_convolution_same_width` uses `gaussian_convolution()` (the `GaussianParams`-based function), not the pdf variant.

### Exponential (inline tests)

**Test:** [`packages/treetime-analytical/src/exponential.rs`](../../packages/treetime-analytical/src/exponential.rs)

**Impl:** [`packages/treetime-analytical/src/exponential.rs`](../../packages/treetime-analytical/src/exponential.rs)

| Test                                          | Purpose                        |
| --------------------------------------------- | ------------------------------ |
| `test_exponential_pdf_positive`               | PDF at positive x              |
| `test_exponential_pdf_negative`               | PDF = 0 for negative x         |
| `test_exponential_pdf_zero`                   | PDF at x=0 equals rate         |
| `test_exponential_convolution_distinct_rates` | Convolution with a != b        |
| `test_exponential_convolution_equal_rates`    | Convolution limit form a = b   |
| `test_exponential_convolution_negative_x`     | Convolution = 0 for negative x |
| `test_exponential_convolution_at_zero`        | Convolution = 0 at x=0         |
| `test_exponential_pdf_grid`                   | PDF grid evaluation            |
| `test_exponential_convolution_grid`           | Convolution grid evaluation    |

### Gaussian-Exponential (inline tests)

**Test:** [`packages/treetime-analytical/src/gaussian_exponential.rs`](../../packages/treetime-analytical/src/gaussian_exponential.rs)

**Impl:** [`packages/treetime-analytical/src/gaussian_exponential.rs`](../../packages/treetime-analytical/src/gaussian_exponential.rs)

| Test                                               | Purpose                      |
| -------------------------------------------------- | ---------------------------- |
| `test_gaussian_exponential_convolution_positive_x` | Convolution at positive x    |
| `test_gaussian_exponential_convolution_grid`       | Grid evaluation non-negative |
| `test_gaussian_exponential_convolution_at_peak`    | Convolution positive at peak |

---

## treetime-ops

**Directory:** [`__tests__/`](../../packages/treetime-ops/src/__tests__/)

**Test:** [`packages/treetime-ops/src/__tests__/multiplication.rs`](../../packages/treetime-ops/src/__tests__/multiplication.rs)

**Impl:** [`packages/treetime-ops/src/multiplication.rs`](../../packages/treetime-ops/src/multiplication.rs)

| Test                                                  | Purpose                       |
| ----------------------------------------------------- | ----------------------------- |
| `test_multiply_many_empty`                            | Empty input returns empty     |
| `test_multiply_many_single`                           | Single array normalization    |
| `test_multiply_many_two_identical`                    | Two identical arrays product  |
| `test_multiply_many_underflow_resistance`             | 100 small values no underflow |
| `test_multiply_many_handles_zero_max`                 | All-zero array handling       |
| `test_lazy_vs_aggressive_same_result_for_few_factors` | Lazy matches aggressive       |
| `test_lazy_normalize_fewer_operations`                | Three arrays lazy normalize   |

### Convolution (inline tests)

**Test:** [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)

**Impl:** [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)

| Test                               | Purpose                         |
| ---------------------------------- | ------------------------------- |
| `test_convolve_riemann_delta`      | Riemann sum with delta function |
| `test_convolve_symmetric`          | Convolution commutativity       |
| `test_convolve_fft_matches_direct` | FFT matches direct convolution  |
