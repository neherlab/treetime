# Supporting Crates Tests

[Back to index](_index.md)

## Summary

| Crate                                       | Module                                                                 | Type                |
| ------------------------------------------- | ---------------------------------------------------------------------- | ------------------- |
| [treetime-primitives](#treetime-primitives) | BitSet128, AsciiChar, Seq                                              | Unit, Parameterized |
| [treetime-utils](#treetime-utils)           | interval, iterator, datetime, fmt, array                               | Unit, Parameterized |
| [treetime-io](#treetime-io)                 | dates_csv, discrete_states_csv, nwk_providers, concat, parse_delimited | Unit, Parameterized |
| [treetime-grid](#treetime-grid)             | Grid, GridFn, interp                                                   | Unit, Parameterized |
| [treetime-cli](#treetime-cli)               | convert (auspice, usher, mutation)                                     | Unit, Parameterized |
| [treetime (alphabet)](#alphabet)            | Alphabet, AlphabetConfig                                               | Unit, Parameterized |
| [treetime (io)](#io)                        | FASTA, Newick                                                          | Unit                |
| [treetime (seq)](#sequence-operations)      | composition, div, char_ranges, mutation                                | Unit, Parameterized |
| [treetime (graph)](#graph)                  | traversal, collapse, edge                                              | Unit                |
| [treetime (repr)](#representation)          | compose_substitutions, discrete_states, marginal_helpers, payloads     | Unit, Parameterized |
| [treetime (prune)](#commands-prune)         | prune, collapse, merge                                                 | Unit                |

---

## treetime-primitives

### BitSet128

**Test:** [`packages/treetime-primitives/src/__tests__/test_bitset128.rs`](../../packages/treetime-primitives/src/__tests__/test_bitset128.rs)

**Impl:** [`packages/treetime-primitives/src/bitset128.rs`](../../packages/treetime-primitives/src/bitset128.rs)

Unit and parameterized tests for all BitSet128 operations.

| Test                                            | Purpose                           |
| ----------------------------------------------- | --------------------------------- |
| `test_bitset128_new`                            | Construct empty set               |
| `test_bitset128_from_iter`                      | Create from char iterator         |
| `test_bitset128_from_slice`                     | Create from char slice            |
| `test_bitset128_is_empty`                       | Empty and non-empty checks        |
| `test_bitset128_len`                            | Element count                     |
| `test_bitset128_clear`                          | Clear all elements                |
| `test_bitset128_insert`                         | Insert with dedup                 |
| `test_bitset128_remove`                         | Remove with idempotence           |
| `test_bitset128_union`                          | Pairwise union                    |
| `test_bitset128_union_with_empty`               | Union with empty set              |
| `test_bitset128_intersection`                   | Pairwise intersection             |
| `test_bitset128_intersection_with_empty`        | Intersection with empty set       |
| `test_bitset128_from_union`                     | Multi-set union                   |
| `test_bitset128_from_union_with_empty`          | Multi-set union with empty        |
| `test_bitset128_from_union_of_all_empty`        | Union of all empty sets           |
| `test_bitset128_from_intersection`              | Multi-set intersection            |
| `test_bitset128_from_intersection_with_empty`   | Multi-set intersection with empty |
| `test_bitset128_from_intersection_of_all_empty` | Intersection of all empty sets    |
| `test_bitset128_difference`                     | Set difference                    |
| `test_bitset128_is_disjoint`                    | Disjoint set check                |
| `test_bitset128_is_subset`                      | Subset relationship check         |
| `test_bitset128_is_superset`                    | Superset relationship check       |
| `test_bitset128_display`                        | Display trait formatting          |
| `test_bitset128_debug`                          | Debug trait formatting            |
| `test_bitset128_eq`                             | Equality and inequality           |
| `test_bitset128_hash`                           | Hash consistency with equality    |
| `test_bitset128_get_empty`                      | Get status for empty set          |
| `test_bitset128_get_unambiguous`                | Get status for single element     |
| `test_bitset128_get_ambiguous`                  | Get status for multiple elements  |
| `test_bitset128_get_one`                        | Get first element from set        |
| `test_bitset128_get_one_exactly`                | Get sole element from singleton   |
| `test_bitset128_add`                            | `+` operator (9 cases)            |
| `test_bitset128_add_ref_right`                  | `+ &rhs` operator (9 cases)       |
| `test_bitset128_add_ref_left`                   | `&lhs +` operator (9 cases)       |
| `test_bitset128_add_ref_both`                   | `&lhs + &rhs` operator (9 cases)  |
| `test_bitset128_add_assign`                     | `+=` operator (9 cases)           |
| `test_bitset128_add_assign_ref`                 | `+= &rhs` operator (9 cases)      |
| `test_bitset128_add_char_to_set`                | `set + u8` operator (4 cases)     |
| `test_bitset128_add_set_to_char`                | `u8 + set` operator (4 cases)     |
| `test_bitset128_or`                             | `\|` operator (9 cases)           |
| `test_bitset128_or_ref_right`                   | `\| &rhs` operator (9 cases)      |
| `test_bitset128_or_ref_left`                    | `&lhs \|` operator (9 cases)      |
| `test_bitset128_or_ref_both`                    | `&lhs \| &rhs` operator (9 cases) |
| `test_bitset128_or_assign`                      | `\|=` operator (9 cases)          |
| `test_bitset128_or_assign_ref`                  | `\|= &rhs` operator (9 cases)     |
| `test_bitset128_or_char_to_set`                 | `set \| u8` operator (4 cases)    |
| `test_bitset128_or_set_to_char`                 | `u8 \| set` operator (4 cases)    |
| `test_bitset128_sub`                            | `-` operator (7 cases)            |
| `test_bitset128_sub_ref_right`                  | `- &rhs` operator (7 cases)       |
| `test_bitset128_sub_ref_left`                   | `&lhs -` operator (7 cases)       |
| `test_bitset128_sub_ref_both`                   | `&lhs - &rhs` operator (7 cases)  |
| `test_bitset128_sub_assign`                     | `-=` operator (7 cases)           |
| `test_bitset128_sub_assign_ref`                 | `-= &rhs` operator (7 cases)      |
| `test_bitset128_sub_char_from_set`              | `set - u8` operator (5 cases)     |
| `test_bitset128_and`                            | `&` operator (10 cases)           |
| `test_bitset128_and_ref_right`                  | `& &rhs` operator (10 cases)      |
| `test_bitset128_and_ref_left`                   | `&lhs &` operator (10 cases)      |
| `test_bitset128_and_ref_both`                   | `&lhs & &rhs` operator (10 cases) |
| `test_bitset128_and_assign`                     | `&=` operator (10 cases)          |
| `test_bitset128_and_assign_ref`                 | `&= &rhs` operator (10 cases)     |
| `test_bitset128_and_char_to_set`                | `set & u8` operator (3 cases)     |
| `test_bitset128_and_set_to_char`                | `u8 & set` operator (3 cases)     |
| `test_bitset128_xor`                            | `^` operator (9 cases)            |
| `test_bitset128_xor_ref_right`                  | `^ &rhs` operator (9 cases)       |
| `test_bitset128_xor_ref_left`                   | `&lhs ^` operator (9 cases)       |
| `test_bitset128_xor_ref_both`                   | `&lhs ^ &rhs` operator (9 cases)  |
| `test_bitset128_xor_assign`                     | `^=` operator (9 cases)           |
| `test_bitset128_xor_assign_ref`                 | `^= &rhs` operator (9 cases)      |
| `test_bitset128_xor_char_to_set`                | `set ^ u8` operator (4 cases)     |
| `test_bitset128_xor_set_to_char`                | `u8 ^ set` operator (4 cases)     |

### AsciiChar

**Test:** [`packages/treetime-primitives/src/__tests__/test_ascii_char.rs`](../../packages/treetime-primitives/src/__tests__/test_ascii_char.rs)

**Impl:** [`packages/treetime-primitives/src/seq_char.rs`](../../packages/treetime-primitives/src/seq_char.rs)

| Test                               | Purpose                          |
| ---------------------------------- | -------------------------------- |
| `test_try_new_valid`               | Construct from valid byte        |
| `test_try_new_invalid_error`       | Reject byte >= 128               |
| `test_try_new_zero_valid`          | Construct from zero byte         |
| `test_try_from_char_valid`         | Construct from ASCII char        |
| `test_try_from_char_invalid_error` | Reject non-ASCII char            |
| `test_deserialize_valid`           | Deserialize valid JSON integer   |
| `test_deserialize_invalid_fails`   | Reject out-of-range JSON integer |

### Seq

**Test:** [`packages/treetime-primitives/src/__tests__/test_seq.rs`](../../packages/treetime-primitives/src/__tests__/test_seq.rs)

**Impl:** [`packages/treetime-primitives/src/seq.rs`](../../packages/treetime-primitives/src/seq.rs)

| Test                                  | Purpose                 |
| ------------------------------------- | ----------------------- |
| `test_try_from_str_non_ascii_error`   | Reject non-ASCII string |
| `test_try_from_vec_non_ascii_error`   | Reject non-ASCII vec    |
| `test_try_from_slice_non_ascii_error` | Reject non-ASCII slice  |

---

## treetime-utils

### interval/range_complement

**Test:** [`packages/treetime-utils/src/interval/range_complement.rs`](../../packages/treetime-utils/src/interval/range_complement.rs) (inline)

**Impl:** [`packages/treetime-utils/src/interval/range_complement.rs`](../../packages/treetime-utils/src/interval/range_complement.rs)

| Test                                           | Purpose                         |
| ---------------------------------------------- | ------------------------------- |
| `test_range_complement_empty`                  | Complement of empty is universe |
| `test_range_complement_multiple`               | Complement of overlapping sets  |
| `test_range_complement_full_universe`          | Full coverage yields empty      |
| `test_range_complement_full_universe_overlap`  | Overlapping full coverage       |
| `test_range_complement_full_universe_adjacent` | Adjacent full coverage          |
| `test_range_complement_full_universe_same`     | Duplicate full coverage         |

### interval/range_difference

**Test:** [`packages/treetime-utils/src/interval/range_difference.rs`](../../packages/treetime-utils/src/interval/range_difference.rs) (inline)

**Impl:** [`packages/treetime-utils/src/interval/range_difference.rs`](../../packages/treetime-utils/src/interval/range_difference.rs)

| Test                                            | Purpose                       |
| ----------------------------------------------- | ----------------------------- |
| `test_range_difference_with_no_overlap`         | No overlap preserves input    |
| `test_range_difference_with_partial_overlap`    | Partial overlap trims range   |
| `test_range_difference_with_full_overlap`       | Full overlap yields empty     |
| `test_range_difference_with_multiple_overlaps`  | Multiple overlaps split range |
| `test_range_difference_with_sub_intervals`      | Sub-intervals punch holes     |
| `test_range_difference_with_adjacent_intervals` | Adjacent interval boundary    |
| `test_range_difference_idempotence`             | A - A = empty                 |
| `test_range_difference_identity`                | A - empty = A                 |

### interval/range_intersection

**Test:** [`packages/treetime-utils/src/interval/range_intersection.rs`](../../packages/treetime-utils/src/interval/range_intersection.rs) (inline)

**Impl:** [`packages/treetime-utils/src/interval/range_intersection.rs`](../../packages/treetime-utils/src/interval/range_intersection.rs)

| Test                                                                 | Purpose                            |
| -------------------------------------------------------------------- | ---------------------------------- |
| `test_range_intersection_empty_input`                                | Empty input yields empty           |
| `test_range_intersection_with_empty`                                 | Intersection with empty set        |
| `test_range_intersection_single_set`                                 | Single set returns itself          |
| `test_range_intersection_multiple_sets_no_overlap`                   | Disjoint sets yield empty          |
| `test_range_intersection_multiple_sets_with_overlap`                 | Overlapping multi-set intersection |
| `test_range_intersection_multiple_sets_with_nested_overlap`          | Nested overlap intersection        |
| `test_range_intersection_multiple_sets_with_same_ranges`             | Identical ranges preserved         |
| `test_range_intersection_single_interval_contained_within_other`     | Contained interval result          |
| `test_range_intersection_single_interval_not_contained_within_other` | Disjoint pair yields empty         |
| `test_range_intersection_multiple_sets_with_partial_overlap`         | Partial overlap narrows result     |
| `test_range_intersection_sets_with_same_start_or_end`                | Shared boundaries preserved        |
| `test_range_intersection_sets_with_same_range_different_count`       | Different set sizes                |
| `test_range_intersection_sets_with_empty_intervals`                  | Empty interval in middle           |
| `test_range_intersection_commutativity`                              | A intersect B = B intersect A      |
| `test_range_intersection_associativity`                              | Associativity property             |
| `test_range_intersection_idempotence`                                | A intersect A = A                  |
| `test_range_intersection_absorption`                                 | Absorption property                |
| `test_range_intersection_empty_set`                                  | Intersection with empty set        |

### interval/range_union

**Test:** [`packages/treetime-utils/src/interval/range_union.rs`](../../packages/treetime-utils/src/interval/range_union.rs) (inline)

**Impl:** [`packages/treetime-utils/src/interval/range_union.rs`](../../packages/treetime-utils/src/interval/range_union.rs)

| Test                                                  | Purpose                       |
| ----------------------------------------------------- | ----------------------------- |
| `test_range_union_empty`                              | Empty input yields empty      |
| `test_range_union_empty_multiple`                     | Multiple empty sets           |
| `test_range_union_with_empty`                         | Union with empty set          |
| `test_range_union_single_set`                         | Single set returns itself     |
| `test_range_union_multiple_sets_no_overlap`           | Disjoint sets concatenated    |
| `test_range_union_multiple_sets_with_overlap`         | Overlapping sets merged       |
| `test_range_union_multiple_sets_with_nested_overlap`  | Nested overlap merged         |
| `test_range_union_multiple_sets_with_same_ranges`     | Identical ranges deduplicated |
| `test_range_union_disjoint_sets`                      | Disjoint multi-range sets     |
| `test_range_union_overlapping_sets_different_lengths` | Different-length overlap      |
| `test_range_union_commutativity`                      | A union B = B union A         |
| `test_range_union_associativity`                      | Associativity property        |
| `test_range_union_idempotence`                        | A union A = A                 |
| `test_range_union_empty_set_identity`                 | A union empty = A             |
| `test_range_union_absorption`                         | Subset absorbed by superset   |

### interval/range_properties

**Test:** [`packages/treetime-utils/src/interval/range_properties.rs`](../../packages/treetime-utils/src/interval/range_properties.rs) (inline)

**Impl:**

- [`packages/treetime-utils/src/interval/range_complement.rs`](../../packages/treetime-utils/src/interval/range_complement.rs)
- [`packages/treetime-utils/src/interval/range_intersection.rs`](../../packages/treetime-utils/src/interval/range_intersection.rs)
- [`packages/treetime-utils/src/interval/range_union.rs`](../../packages/treetime-utils/src/interval/range_union.rs)

| Test                                                         | Purpose                                   |
| ------------------------------------------------------------ | ----------------------------------------- |
| `test_range_properties_union_distribution_over_intersection` | Distributivity of union over intersection |
| `test_range_properties_intersection_absorbs_union`           | Intersection absorbs union                |
| `test_range_properties_union_absorbs_intersection`           | Union absorbs intersection                |
| `test_range_union_with_complement`                           | A union complement(A) = universe          |
| `test_range_union_with_complement_multiple_sets`             | Multi-set complement union                |
| `test_range_union_complement_empty_set`                      | Complement of empty = universe            |

### iterator/mean_by_key

**Test:** [`packages/treetime-utils/src/iterator/mean_by_key.rs`](../../packages/treetime-utils/src/iterator/mean_by_key.rs) (inline)

**Impl:** [`packages/treetime-utils/src/iterator/mean_by_key.rs`](../../packages/treetime-utils/src/iterator/mean_by_key.rs)

| Test                                      | Purpose                     |
| ----------------------------------------- | --------------------------- |
| `test_mean_by_key_i32`                    | Mean with i32 transform     |
| `test_mean_by_key_f64`                    | Mean of f64 values          |
| `test_mean_by_key_f32`                    | Mean with f32 transform     |
| `test_mean_by_key_empty`                  | Empty iterator returns zero |
| `test_mean_by_key_single`                 | Single-element mean         |
| `test_mean_by_key_struct`                 | Mean of struct fields       |
| `test_mean_by_key_owned`                  | Owned iterator mean         |
| `test_mean_by_key_negative_numbers`       | Negative integers cancel    |
| `test_mean_by_key_negative_floats`        | Negative floats cancel      |
| `test_mean_by_key_large_collection`       | 1000-element mean           |
| `test_mean_by_key_all_zeros`              | All-zero input              |
| `test_mean_by_key_u32`                    | Mean of u32 values          |
| `test_mean_by_key_i64`                    | Mean of large i64 values    |
| `test_mean_by_key_complex_transformation` | Quadratic transform mean    |
| `test_mean_by_key_very_small_floats`      | Near-zero float mean        |
| `test_mean_by_key_mixed_precision`        | Decimal precision mean      |
| `test_mean_by_key_chained_iterator`       | Filtered iterator mean      |
| `test_mean_by_key_absolute_values`        | Absolute value transform    |
| `test_mean_by_key_option_values`          | Filter-mapped Option mean   |

### datetime/parse_uncertain_date

**Test:** [`packages/treetime-utils/src/datetime/__tests__/parse_uncertain_date.rs`](../../packages/treetime-utils/src/datetime/__tests__/parse_uncertain_date.rs)

**Impl:** [`packages/treetime-utils/src/datetime/parse_uncertain_date.rs`](../../packages/treetime-utils/src/datetime/parse_uncertain_date.rs)

| Test                              | Purpose                                         |
| --------------------------------- | ----------------------------------------------- |
| `test_date_parse_uncertain`       | Parse date with wildcards (41 cases)            |
| `test_date_parse_uncertain_error` | Reject fully-unknown and empty dates (14 cases) |

### datetime/format_to_regex

**Test:** [`packages/treetime-utils/src/datetime/format_to_regex.rs`](../../packages/treetime-utils/src/datetime/format_to_regex.rs) (inline)

**Impl:** [`packages/treetime-utils/src/datetime/format_to_regex.rs`](../../packages/treetime-utils/src/datetime/format_to_regex.rs)

| Test                                  | Purpose                             |
| ------------------------------------- | ----------------------------------- |
| `test_date_format_to_regex`           | Convert strftime to regex (7 cases) |
| `test_date_format_to_regex_uncertain` | Match uncertain XX dates (2 cases)  |

### datetime/year_fraction

**Test:** [`packages/treetime-utils/src/datetime/year_fraction.rs`](../../packages/treetime-utils/src/datetime/year_fraction.rs) (inline)

**Impl:** [`packages/treetime-utils/src/datetime/year_fraction.rs`](../../packages/treetime-utils/src/datetime/year_fraction.rs)

| Test                         | Purpose                                        |
| ---------------------------- | ---------------------------------------------- |
| `test_date_to_year_fraction` | DateTime to year fraction roundtrip (50 cases) |

### fmt/string

**Test:** [`packages/treetime-utils/src/fmt/__tests__/string.rs`](../../packages/treetime-utils/src/fmt/__tests__/string.rs)

**Impl:** [`packages/treetime-utils/src/fmt/string.rs`](../../packages/treetime-utils/src/fmt/string.rs)

| Test                                 | Purpose                                  |
| ------------------------------------ | ---------------------------------------- |
| `test_truncate_right`                | Truncate from right (16 cases)           |
| `test_truncate_left`                 | Truncate from left (16 cases)            |
| `test_truncate_middle`               | Truncate from middle (16 cases)          |
| `test_truncate_right_with_ellipsis`  | Right truncate with `...` (16 cases)     |
| `test_truncate_left_with_ellipsis`   | Left truncate with `...` (19 cases)      |
| `test_truncate_middle_with_ellipsis` | Middle truncate with `...` (19 cases)    |
| `test_truncate_with_custom_ellipsis` | Custom ellipsis and direction (16 cases) |

### fmt/float

**Test:** [`packages/treetime-utils/src/fmt/float.rs`](../../packages/treetime-utils/src/fmt/float.rs) (inline)

**Impl:** [`packages/treetime-utils/src/fmt/float.rs`](../../packages/treetime-utils/src/fmt/float.rs)

| Test                               | Purpose                                   |
| ---------------------------------- | ----------------------------------------- |
| `test_float_to_significant_digits` | Format to significant digits (4 cases)    |
| `test_float_to_decimal_digits`     | Format to decimal digits (4 cases)        |
| `test_float_to_digits`             | Format with both constraints (4 cases)    |
| `test_trait_to_significant_digits` | Trait method significant digits (3 cases) |
| `test_trait_to_decimal_digits`     | Trait method decimal digits (3 cases)     |
| `test_trait_to_digits`             | Trait method both constraints (4 cases)   |
| `test_trim_trailing_zeros`         | Trim trailing zeros (8 cases)             |
| `test_edge_cases`                  | Zero, negative, infinity (4 cases)        |

### array/ndarray

**Test:** [`packages/treetime-utils/src/array/__tests__/ndarray.rs`](../../packages/treetime-utils/src/array/__tests__/ndarray.rs)

**Impl:** [`packages/treetime-utils/src/array/ndarray.rs`](../../packages/treetime-utils/src/array/ndarray.rs)

| Test                                         | Purpose                           |
| -------------------------------------------- | --------------------------------- |
| `computes_argmin_axis_0`                     | Argmin along axis 0               |
| `computes_argmin_axis_1`                     | Argmin along axis 1               |
| `computes_argmax_axis_0`                     | Argmax along axis 0               |
| `computes_argmax_axis_1`                     | Argmax along axis 1               |
| `computes_cumsum_axis_0`                     | Cumulative sum along axis 0       |
| `computes_cumsum_axis_1`                     | Cumulative sum along axis 1       |
| `computes_outer_product`                     | Basic outer product               |
| `test_outer_product_case_1`                  | Outer product with uniform pi     |
| `test_outer_product_case_2`                  | Outer product with non-uniform pi |
| `chooses_1d_by_index`                        | Index-based 1D selection          |
| `chooses_2d_by_index`                        | Index-based 2D selection          |
| `generates_predictable_random_uniform`       | Seeded random matrix              |
| `test_product_axis_empty_array1`             | Product of empty 1D array         |
| `test_product_axis_empty_array2_axis0`       | Product of empty 2D along axis 0  |
| `test_product_axis_empty_array2_axis1`       | Product of empty 2D along axis 1  |
| `test_product_axis_empty_axis_0_prod_axis_0` | Single-row product axis 0         |
| `test_product_axis_empty_axis_0_prod_axis_1` | Single-row product axis 1         |
| `test_product_axis_empty_axis_1_prod_axis_0` | Single-column product axis 0      |
| `test_product_axis_empty_axis_1_prod_axis_1` | Single-column product axis 1      |
| `test_product_axis_general_case_axis_0`      | General product axis 0            |
| `test_product_axis_general_case_axis_1`      | General product axis 1            |

---

## treetime-io

### dates_csv

**Test:** [`packages/treetime-io/src/__tests__/test_dates_csv.rs`](../../packages/treetime-io/src/__tests__/test_dates_csv.rs)

**Impl:** [`packages/treetime-io/src/dates_csv.rs`](../../packages/treetime-io/src/dates_csv.rs)

| Test                       | Purpose                                  |
| -------------------------- | ---------------------------------------- |
| `test_date_read_empty`     | Parse empty/null date strings (7 cases)  |
| `test_date_read_ok`        | Parse dates to year fractions (12 cases) |
| `test_read_dates_from_str` | Parse TSV metadata file                  |

### discrete_states_csv

**Test:** [`packages/treetime-io/src/__tests__/test_discrete_states_csv.rs`](../../packages/treetime-io/src/__tests__/test_discrete_states_csv.rs)

**Impl:** [`packages/treetime-io/src/discrete_states_csv.rs`](../../packages/treetime-io/src/discrete_states_csv.rs)

| Test                                                        | Purpose                                               |
| ----------------------------------------------------------- | ----------------------------------------------------- |
| `test_discrete_states_csv_delimiter_handling` (2 cases)     | TSV and CSV delimiters parse correctly                |
| `test_discrete_states_csv_name_column_detection` (4 cases)  | Default and custom name column headers detected       |
| `test_discrete_states_csv_value_column_selection` (2 cases) | Explicit value column selected among multiple columns |
| `test_discrete_states_csv_header_normalization`             | Headers with `#` delimiters normalized correctly      |
| `test_discrete_states_csv_custom_parser`                    | Custom value parser (f64 parsing) works               |
| `test_discrete_states_csv_whitespace_trimming`              | Whitespace in names and values trimmed                |

### nwk_providers

**Test:** [`packages/treetime-io/src/__tests__/test_nwk_providers.rs`](../../packages/treetime-io/src/__tests__/test_nwk_providers.rs)

**Impl:** [`packages/treetime-io/src/nwk.rs`](../../packages/treetime-io/src/nwk.rs)

| Test                                        | Purpose                                                                    |
| ------------------------------------------- | -------------------------------------------------------------------------- |
| `test_comment_provider_empty`               | Empty CommentProviders produces same output as default nwk_write_str       |
| `test_comment_provider_single`              | Single provider adds NHX-style comment to one node                         |
| `test_comment_provider_multiple`            | Multiple providers append separate comment blocks                          |
| `test_comment_provider_precedence`          | Later provider overwrites earlier provider for same key                    |
| `test_comment_provider_overwrites_payload`  | Provider comments overwrite payload comments for same key, preserve others |
| `test_comment_provider_serializes_comments` | Provider comments serialize into Newick comment format                     |

### concat

**Test:** [`packages/treetime-io/src/__tests__/test_concat.rs`](../../packages/treetime-io/src/__tests__/test_concat.rs)

**Impl:** [`packages/treetime-io/src/concat.rs`](../../packages/treetime-io/src/concat.rs)

| Test                                                              | Purpose                                |
| ----------------------------------------------------------------- | -------------------------------------- |
| `test_concatenate_with_delimiter_both_with_trailing_newline`      | Two readers with trailing newlines     |
| `test_concatenate_with_delimiter_one_with_trailing_newline`       | Single reader with trailing newline    |
| `test_concatenate_with_delimiter_one_without_trailing_newline`    | Single reader without trailing newline |
| `test_concatenate_with_delimiter_one_empty`                       | Single empty reader                    |
| `test_concatenate_with_delimiter_one_empty_with_trailing_newline` | Single newline-only reader             |
| `test_concatenate_with_delimiter_no_newline`                      | Two readers without trailing newlines  |
| `test_concatenate_with_delimiter_first_empty_no_newline`          | First reader empty                     |
| `test_concatenate_with_delimiter_second_empty_no_newline`         | Second reader empty                    |
| `test_concatenate_with_delimiter_both_empty`                      | Both readers empty                     |

### parse_delimited

**Test:** [`packages/treetime-io/src/parse_delimited.rs`](../../packages/treetime-io/src/parse_delimited.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime-io/src/parse_delimited.rs`](../../packages/treetime-io/src/parse_delimited.rs)

| Test                                                         | Purpose                                          |
| ------------------------------------------------------------ | ------------------------------------------------ |
| `test_parse_delimited_str_comma_delimiter`                   | Comma-separated values parsed                    |
| `test_parse_delimited_str_semicolon_delimiter`               | Semicolon-separated values parsed                |
| `test_parse_delimited_str_pipe_delimiter`                    | Pipe-separated values parsed                     |
| `test_parse_delimited_str_newline_delimiter`                 | Newline-separated values parsed                  |
| `test_parse_delimited_str_whitespace_in_names`               | Whitespace preserved within elements             |
| `test_parse_delimited_str_empty_string_produces_empty_vec`   | Empty input produces empty result                |
| `test_parse_delimited_str_single_element_no_delimiter`       | Single element without delimiter parsed          |
| `test_parse_delimited_str_empty_elements_between_delimiters` | Empty elements between consecutive delimiters    |
| `test_parse_delimited_str_trailing_delimiter`                | Trailing delimiter behavior                      |
| `test_parse_delimited_str_leading_delimiter`                 | Leading delimiter produces leading empty element |
| `test_parse_delimited_str_only_delimiters`                   | Only-delimiter input behavior                    |
| `test_parse_delimited_str_tab_delimiter`                     | Tab-separated values parsed                      |
| `test_parse_delimited_str_space_delimiter`                   | Space-separated values parsed                    |
| `test_parse_delimited_str_unicode_content`                   | Unicode content preserved                        |
| `test_parse_delimited_str_mixed_whitespace_content`          | Mixed whitespace content preserved               |

---

## treetime-grid

### Grid

**Test:** [`packages/treetime-grid/src/__tests__/grid.rs`](../../packages/treetime-grid/src/__tests__/grid.rs)

**Impl:** [`packages/treetime-grid/src/grid.rs`](../../packages/treetime-grid/src/grid.rs)

| Test                       | Purpose                                       |
| -------------------------- | --------------------------------------------- |
| `test_find_interval_index` | Find grid interval for query point (11 cases) |

### GridFn

**Test:** [`packages/treetime-grid/src/__tests__/grid_fn.rs`](../../packages/treetime-grid/src/__tests__/grid_fn.rs)

**Impl:** [`packages/treetime-grid/src/grid_fn.rs`](../../packages/treetime-grid/src/grid_fn.rs)

| Test                                     | Purpose                                    |
| ---------------------------------------- | ------------------------------------------ |
| `test_gridfn_interp_exact_grid_points`   | Interpolate at exact grid points (5 cases) |
| `test_gridfn_interp_interior_points`     | Interpolate between grid points (7 cases)  |
| `test_gridfn_interp_left_extrapolation`  | Constant extrapolation below min (5 cases) |
| `test_gridfn_interp_right_extrapolation` | Constant extrapolation above max (5 cases) |
| `test_gridfn_interp_many`                | Batch interpolation                        |
| `test_gridfn_from_grid`                  | Construct from function                    |
| `test_gridfn_accessors`                  | X and Y accessor methods                   |
| `test_gridfn_mapv`                       | Map values transform                       |
| `test_gridfn_negate_arg_inplace`         | Negate argument in-place                   |
| `test_gridfn_resample_to_grid_finer`     | Resample to finer grid                     |

### interp_nonuniform

**Test:** [`packages/treetime-grid/src/interp_nonuniform.rs`](../../packages/treetime-grid/src/interp_nonuniform.rs) (inline)

**Impl:** [`packages/treetime-grid/src/interp_nonuniform.rs`](../../packages/treetime-grid/src/interp_nonuniform.rs)

| Test                                         | Purpose                                |
| -------------------------------------------- | -------------------------------------- |
| `test_interp_nonuniform_exact_points`        | Interpolate at known points            |
| `test_interp_nonuniform_interior`            | Interpolate between non-uniform points |
| `test_interp_nonuniform_left_extrapolation`  | Linear extrapolation below min         |
| `test_interp_nonuniform_right_extrapolation` | Linear extrapolation above max         |
| `test_interp_nonuniform_mixed`               | Mixed interpolation and extrapolation  |

---

## treetime-cli

### convert/auspice

**Test:** [`packages/treetime-cli/src/convert/__tests__/auspice.rs`](../../packages/treetime-cli/src/convert/__tests__/auspice.rs)

**Impl:** [`packages/treetime-cli/src/convert/auspice.rs`](../../packages/treetime-cli/src/convert/auspice.rs)

| Test                                      | Purpose                            |
| ----------------------------------------- | ---------------------------------- |
| `test_auspice_from_nwk`                   | Convert Newick to Auspice JSON     |
| `test_auspice_to_nwk`                     | Convert Auspice JSON to Newick     |
| `test_auspice_roundtrip`                  | Auspice JSON read-write roundtrip  |
| `test_auspice_mutation_roundtrip`         | Mutation preservation in roundtrip |
| `test_nwk_auspice_nwk_topology_preserved` | Topology preserved across formats  |

### convert/convert

**Test:** [`packages/treetime-cli/src/convert/__tests__/convert.rs`](../../packages/treetime-cli/src/convert/__tests__/convert.rs)

**Impl:**

- [`packages/treetime-cli/src/convert/convert.rs`](../../packages/treetime-cli/src/convert/convert.rs)
- [`packages/treetime-cli/src/convert/auspice.rs`](../../packages/treetime-cli/src/convert/auspice.rs)
- [`packages/treetime-cli/src/convert/usher.rs`](../../packages/treetime-cli/src/convert/usher.rs)

| Test                                         | Purpose                            |
| -------------------------------------------- | ---------------------------------- |
| `test_newick_to_auspice_to_newick_topology`  | Newick-Auspice-Newick roundtrip    |
| `test_auspice_to_newick_to_auspice_topology` | Auspice-Newick-Auspice roundtrip   |
| `test_auspice_to_usher_mutation_transfer`    | Auspice to UShER mutation transfer |

### convert/usher

**Test:** [`packages/treetime-cli/src/convert/__tests__/usher.rs`](../../packages/treetime-cli/src/convert/__tests__/usher.rs)

**Impl:** [`packages/treetime-cli/src/convert/usher.rs`](../../packages/treetime-cli/src/convert/usher.rs)

| Test                                       | Purpose                            |
| ------------------------------------------ | ---------------------------------- |
| `test_usher_write_mutations`               | Write nuc mutations to UShER       |
| `test_auspice_usher_auspice_nuc_mutations` | Only nuc mutations in UShER output |

### convert/mutation

**Test:** [`packages/treetime-cli/src/convert/mutation.rs`](../../packages/treetime-cli/src/convert/mutation.rs) (inline)

**Impl:** [`packages/treetime-cli/src/convert/mutation.rs`](../../packages/treetime-cli/src/convert/mutation.rs)

| Test                        | Purpose                                   |
| --------------------------- | ----------------------------------------- |
| `test_mutation_parse`       | Parse mutation strings (5 cases)          |
| `test_mutation_parse_error` | Reject invalid mutation strings (6 cases) |
| `test_mutation_format`      | Format mutation to string (3 cases)       |
| `test_mutation_roundtrip`   | Parse-format roundtrip (4 cases)          |
| `test_parse_mutation_list`  | Parse list of mutations                   |
| `test_format_mutation_list` | Format list of mutations                  |

---

## treetime Core Library

### Alphabet

**Test:** [`packages/treetime/src/alphabet/__tests__/test_alphabet.rs`](../../packages/treetime/src/alphabet/__tests__/test_alphabet.rs)

**Impl:**

- [`packages/treetime/src/alphabet/alphabet.rs`](../../packages/treetime/src/alphabet/alphabet.rs)
- [`packages/treetime/src/alphabet/alphabet_config.rs`](../../packages/treetime/src/alphabet/alphabet_config.rs)

| Test                                             | Purpose                                      |
| ------------------------------------------------ | -------------------------------------------- |
| `test_alphabet_default`                          | Default alphabet has 4 canonical             |
| `test_alphabet_new_n_canonical` (2 cases)        | Nuc has 4, Aa has 21 canonical chars         |
| `test_alphabet_new_n_ambiguous` (2 cases)        | Nuc has 10, Aa has 3 ambiguous chars         |
| `test_alphabet_new_n_undetermined` (2 cases)     | Nuc has 2, Aa has 2 undetermined chars       |
| `test_alphabet_unknown` (2 cases)                | Nuc unknown is N, Aa unknown is X            |
| `test_alphabet_gap` (2 cases)                    | Nuc/Aa gap char is dash                      |
| `test_alphabet_nuc_canonical_chars`              | Nuc canonical is A,C,G,T                     |
| `test_alphabet_aa_canonical_chars`               | Aa canonical is 21 amino acids               |
| `test_alphabet_nuc_is_canonical` (6 cases)       | A,C,G,T canonical; N,gap not                 |
| `test_alphabet_nuc_is_ambiguous` (7 cases)       | R,Y,S,W ambiguous; A,N,gap not               |
| `test_alphabet_nuc_is_determined` (4 cases)      | A,R determined; N,gap not                    |
| `test_alphabet_nuc_is_undetermined` (4 cases)    | N,gap undetermined; A,R not                  |
| `test_alphabet_nuc_is_unknown` (3 cases)         | N unknown; gap,A not                         |
| `test_alphabet_nuc_is_gap` (3 cases)             | Gap is gap; N,A not                          |
| `test_alphabet_nuc_contains` (6 cases)           | A,N,gap,R in; X,Z not                        |
| `test_alphabet_nuc_n_chars`                      | Nuc alphabet has 16 total chars              |
| `test_alphabet_nuc_n_determined`                 | Nuc has 14 determined chars                  |
| `test_alphabet_nuc_index` (4 cases)              | A->0, C->1, G->2, T->3                       |
| `test_alphabet_nuc_char` (4 cases)               | 0->A, 1->C, 2->G, 3->T                       |
| `test_alphabet_nuc_get_profile_canonical`        | Canonical chars produce identity profiles    |
| `test_alphabet_nuc_get_profile_ambiguous`        | Ambiguous chars produce multi-state profiles |
| `test_alphabet_nuc_get_profile_unknown`          | Unknown char produces all-ones profile       |
| `test_alphabet_nuc_get_code`                     | Profile-to-char reverse lookup               |
| `test_alphabet_nuc_char_to_set_canonical`        | Canonical char maps to singleton set         |
| `test_alphabet_nuc_char_to_set_ambiguous`        | R maps to set of A and G                     |
| `test_alphabet_nuc_set_to_char`                  | Set-to-char reverse lookup                   |
| `test_alphabet_nuc_construct_profile`            | Profile from A,G pair                        |
| `test_alphabet_nuc_seq2prof`                     | Sequence to profile matrix                   |
| `test_alphabet_gap_is_not_unknown`               | Gap distinct from unknown                    |
| `test_alphabet_gap_profile_matches_unknown`      | Gap profile matches unknown profile          |
| `test_alphabet_reserved_constants`               | NON_CHAR, VARIABLE_CHAR, FILL_CHAR values    |
| `test_alphabet_aa_ambiguous`                     | Aa ambiguous chars are B,J,Z                 |
| `test_alphabet_aa_b_maps_to_nd`                  | Aa B maps to N,D set                         |
| `test_alphabet_with_config_empty_canonical`      | Empty canonical config errors                |
| `test_alphabet_nuc_chars_iterator`               | All-chars iterator yields 16                 |
| `test_alphabet_nuc_determined_iterator`          | Determined iterator yields 14                |
| `test_alphabet_nuc_undetermined_iterator`        | Undetermined iterator yields N,gap           |
| `test_alphabet_nuc_ambiguous_iterator`           | Ambiguous iterator yields 10                 |
| `test_alphabet_construct_profile_with_ambiguous` | Profile from ambiguous R input               |
| `test_alphabet_construct_profile_empty`          | Empty input produces zero profile            |
| `test_alphabet_seq2prof_with_ambiguous`          | Seq2prof handles ambiguous chars             |
| `test_alphabet_name_display` (2 cases)           | Nuc display is "Nuc", Aa is "Aa"             |
| `test_alphabet_name_default`                     | Default alphabet name is Nuc                 |
| `test_alphabet_with_custom_config`               | Custom XYZ alphabet construction             |
| `test_alphabet_nuc_all_ambiguous_mappings`       | All 10 nuc IUPAC ambiguity mappings          |
| `test_alphabet_aa_all_ambiguous_mappings`        | All 3 aa ambiguity mappings                  |
| `test_alphabet_char_index_roundtrip`             | Char-to-index-to-char roundtrip              |
| `test_alphabet_set_to_char_canonical_roundtrip`  | Char-to-set-to-char roundtrip                |

**Test:** [`packages/treetime/src/alphabet/__tests__/test_alphabet_config.rs`](../../packages/treetime/src/alphabet/__tests__/test_alphabet_config.rs)

**Impl:** [`packages/treetime/src/alphabet/alphabet_config.rs`](../../packages/treetime/src/alphabet/alphabet_config.rs)

| Test                                                                    | Purpose                                           |
| ----------------------------------------------------------------------- | ------------------------------------------------- |
| `test_alphabet_config_validate_valid`                                   | Valid config passes validation                    |
| `test_alphabet_config_validate_duplicate_canonical`                     | Duplicate canonical chars rejected                |
| `test_alphabet_config_validate_canonical_ambiguous_overlap`             | Canonical-ambiguous overlap rejected              |
| `test_alphabet_config_validate_canonical_contains_gap`                  | Gap in canonical rejected                         |
| `test_alphabet_config_validate_canonical_contains_unknown`              | Unknown in canonical rejected                     |
| `test_alphabet_config_validate_ambiguous_contains_gap`                  | Gap as ambiguous key rejected                     |
| `test_alphabet_config_validate_ambiguous_maps_to_noncanonical`          | Ambiguous mapping to non-canonical rejected       |
| `test_alphabet_config_validate_reserved_in_canonical` (3 cases)         | Reserved chars in canonical rejected              |
| `test_alphabet_config_validate_reserved_in_ambiguous_key` (3 cases)     | Reserved chars as ambiguous key rejected          |
| `test_alphabet_config_validate_reserved_in_ambiguous_value` (3 cases)   | Reserved chars in ambiguous value rejected        |
| `test_alphabet_config_create_profile_map`                               | Profile map for canonical and ambiguous           |
| `test_alphabet_config_create_profile_map_gap_matches_unknown`           | Gap profile matches unknown profile               |
| `test_alphabet_config_create_profile_map_identity_matrix_for_canonical` | Canonical chars produce identity matrix           |
| `test_alphabet_config_serde_roundtrip`                                  | JSON serialize-deserialize roundtrip              |
| `test_alphabet_config_validate_reserved_as_unknown` (3 cases)           | Reserved chars as unknown rejected                |
| `test_alphabet_config_validate_reserved_as_gap` (3 cases)               | Reserved chars as gap rejected                    |
| `test_alphabet_config_validate_ambiguous_empty_values`                  | Empty ambiguous values pass validation            |
| `test_alphabet_config_validate_empty_canonical`                         | Empty canonical passes config validation          |
| `test_alphabet_config_create_profile_map_all_ambiguous`                 | Multiple ambiguous codes produce correct profiles |
| `test_alphabet_config_equality`                                         | Config equality and inequality                    |
| `test_alphabet_config_clone`                                            | Config clone matches original                     |
| `test_alphabet_config_validate_gap_equals_unknown`                      | Gap equal to unknown passes validation            |

### I/O

**Test:** [`packages/treetime/src/io/__tests__/test_fasta.rs`](../../packages/treetime/src/io/__tests__/test_fasta.rs)

**Impl:** [`packages/treetime-io/src/fasta.rs`](../../packages/treetime-io/src/fasta.rs)

| Test                                                                  | Purpose                                       |
| --------------------------------------------------------------------- | --------------------------------------------- |
| `test_fasta_reader_fail_on_non_fasta`                                 | Non-FASTA input errors                        |
| `test_fasta_reader_fail_on_unknown_char`                              | Unknown character in sequence errors          |
| `test_fasta_reader_read_empty`                                        | Empty input yields empty record               |
| `test_fasta_reader_read_whitespace_only`                              | Whitespace-only input yields empty record     |
| `test_fasta_reader_read_single_record`                                | Single record parsed correctly                |
| `test_fasta_reader_read_single_record_with_leading_newline`           | Leading newline before first record           |
| `test_fasta_reader_read_single_record_with_multiple_leading_newlines` | Multiple leading newlines handled             |
| `test_fasta_reader_read_single_record_without_trailing_newline`       | Missing trailing newline handled              |
| `test_fasta_reader_read_multiple_records`                             | Two records parsed sequentially               |
| `test_fasta_reader_read_empty_lines_between_records`                  | Empty lines between records handled           |
| `test_fasta_reader_read_with_trailing_newline`                        | Trailing newline after sequence               |
| `test_fasta_reader_example_1`                                         | Two records without trailing newline          |
| `test_fasta_reader_example_2`                                         | Three records including empty sequence        |
| `test_fasta_reader_example_3`                                         | Back-to-back headers with empty sequence      |
| `test_fasta_reader_name_desc`                                         | Name and description parsing                  |
| `test_fasta_reader_dedent_nuc`                                        | Multiple nuc records with varied formatting   |
| `test_fasta_reader_dedent_aa`                                         | Multiple amino acid records                   |
| `test_fasta_reader_multiline_and_skewed_indentation`                  | Case folding, multiline, indentation handling |

**Test:** [`packages/treetime/src/io/__tests__/test_nwk.rs`](../../packages/treetime/src/io/__tests__/test_nwk.rs)

**Impl:** [`packages/treetime-io/src/nwk.rs`](../../packages/treetime-io/src/nwk.rs)

| Test                                      | Purpose                              |
| ----------------------------------------- | ------------------------------------ |
| `test_nwk_roundtrip_binary_tree`          | Binary tree Newick roundtrip         |
| `test_nwk_parse_no_branch_lengths`        | Missing branch lengths parsed as NaN |
| `test_nwk_roundtrip_single_node`          | Single node tree roundtrip           |
| `test_nwk_roundtrip_polytomy`             | Three-child polytomy roundtrip       |
| `test_nwk_roundtrip_nested_polytomy`      | Nested polytomies roundtrip          |
| `test_nwk_roundtrip_zero_length_branches` | Zero-length branches preserved       |
| `test_nwk_parse_verifies_branch_lengths`  | Branch length values match input     |
| `test_nwk_parse_verifies_leaf_names`      | Leaf names match input               |

### Sequence Operations

**Test:** [`packages/treetime/src/seq/__tests__/test_composition.rs`](../../packages/treetime/src/seq/__tests__/test_composition.rs)

**Impl:** [`packages/treetime/src/seq/composition.rs`](../../packages/treetime/src/seq/composition.rs)

| Test                                          | Purpose                              |
| --------------------------------------------- | ------------------------------------ |
| `test_composition_empty`                      | Empty composition from char set      |
| `test_composition_with_sequence`              | Composition from nucleotide sequence |
| `test_composition_with_sequence_and_alphabet` | Composition with full alphabet       |
| `test_composition_add_sequence`               | Incremental sequence addition        |
| `test_composition_add_sequence_with_refs`     | Add sequence from Seq reference      |
| `test_composition_add_mutation`               | Substitution updates counts          |
| `test_composition_add_deletion`               | Deletion updates gap counts          |
| `test_composition_add_insertion`              | Insertion updates char counts        |

**Test:** [`packages/treetime/src/seq/__tests__/test_div.rs`](../../packages/treetime/src/seq/__tests__/test_div.rs)

**Impl:** [`packages/treetime/src/seq/div.rs`](../../packages/treetime/src/seq/div.rs)

| Test                       | Purpose                              |
| -------------------------- | ------------------------------------ |
| `test_all_nodes`           | Root-to-tip divergence for all nodes |
| `test_only_leaves`         | Root-to-tip divergence leaves only   |
| `test_unnamed_internals`   | Unnamed internal nodes handled       |
| `test_single_node`         | Single node has zero divergence      |
| `test_linear_chain`        | Linear chain divergence sum          |
| `test_deep_tree`           | 20-level deep tree divergence        |
| `test_zero_branch_lengths` | Zero branch lengths handled          |

**Test:** [`packages/treetime/src/seq/__tests__/test_find_char_ranges.rs`](../../packages/treetime/src/seq/__tests__/test_find_char_ranges.rs)

**Impl:** [`packages/treetime/src/seq/find_char_ranges.rs`](../../packages/treetime/src/seq/find_char_ranges.rs)

| Test                                       | Purpose                                  |
| ------------------------------------------ | ---------------------------------------- |
| `test_find_letter_ranges` (7 cases)        | Letter range detection in sequences      |
| `test_find_ambiguous_ranges` (7 cases)     | Ambiguous (N) range detection            |
| `test_find_gap_ranges` (7 cases)           | Gap range detection                      |
| `test_find_undetermined_ranges` (12 cases) | Combined N+gap range detection and merge |

**Test:** [`packages/treetime/src/seq/__tests__/test_mutation.rs`](../../packages/treetime/src/seq/__tests__/test_mutation.rs)

**Impl:** [`packages/treetime/src/seq/mutation.rs`](../../packages/treetime/src/seq/mutation.rs)

| Test                                                            | Purpose                                          |
| --------------------------------------------------------------- | ------------------------------------------------ |
| `test_mutation_compose_substitutions_both_empty`                | Empty parent and child produce empty result      |
| `test_mutation_compose_substitutions_parent_empty`              | Empty parent passes child through                |
| `test_mutation_compose_substitutions_child_empty`               | Empty child passes parent through                |
| `test_mutation_compose_substitutions_non_overlapping`           | Disjoint positions merge sorted                  |
| `test_mutation_compose_substitutions_chain`                     | A->G + G->T = A->T at same position              |
| `test_mutation_compose_substitutions_cancellation`              | A->G + G->A = no mutation                        |
| `test_mutation_compose_substitutions_mixed`                     | Chain, passthrough, and cancellation in one call |
| `test_mutation_compose_substitutions_output_sorted_by_position` | Interleaved positions verify merge order         |
| `test_mutation_compose_substitutions_all_cancel`                | All positions cancel, empty result               |

### Graph

**Test:** [`packages/treetime/src/graph/__tests__/graph.rs`](../../packages/treetime/src/graph/__tests__/graph.rs)

**Impl:**

- [`packages/treetime-graph/src/graph.rs`](../../packages/treetime-graph/src/graph.rs)
- [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs)
- [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)

| Test                                                  | Purpose                                |
| ----------------------------------------------------- | -------------------------------------- |
| `test_traversal_serial_depth_first_preorder_forward`  | DFS preorder visits root first         |
| `test_traversal_serial_depth_first_postorder_forward` | DFS postorder visits leaves first      |
| `test_traversal_serial_breadth_first_forward`         | BFS forward visits by level            |
| `test_traversal_serial_breadth_first_reverse`         | BFS reverse visits deepest first       |
| `test_traversal_parallel_breadth_first_forward`       | Parallel BFS forward order             |
| `test_traversal_parallel_breadth_first_backward`      | Parallel BFS backward order            |
| `test_collapse_edge_simple_chain`                     | Collapse single edge in chain          |
| `test_collapse_edge_binary_tree`                      | Collapse internal node in binary tree  |
| `test_collapse_edge_complex_tree`                     | Collapse promotes children to parent   |
| `test_collapse_edge_invalid_edge`                     | Invalid edge key errors                |
| `test_collapse_edge_leaf_edge`                        | Collapse edge to leaf removes leaf     |
| `test_collapse_edge_no_duplicate_edges`               | No duplicate edges after collapse      |
| `test_collapse_edge_adjacency_lists_maintained`       | Adjacency lists correct after collapse |
| `test_collapse_edge_multiple_inbound_edges`           | Collapse with multiple inbound edges   |
| `test_collapse_edge_adjacency_consistency`            | Edge-node adjacency stays consistent   |

**Test:** [`packages/treetime/src/graph/__tests__/test_edge.rs`](../../packages/treetime/src/graph/__tests__/test_edge.rs)

**Impl:** [`packages/treetime-graph/src/edge.rs`](../../packages/treetime-graph/src/edge.rs)

| Test           | Purpose                                |
| -------------- | -------------------------------------- |
| `edge_inverts` | Edge inversion swaps source and target |

### Representation

**Test:** [`packages/treetime/src/representation/__tests__/test_partition_marginal_sparse.rs`](../../packages/treetime/src/representation/__tests__/test_partition_marginal_sparse.rs)

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

**Test:** [`packages/treetime/src/representation/discrete_states.rs`](../../packages/treetime/src/representation/discrete_states.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/representation/discrete_states.rs`](../../packages/treetime/src/representation/discrete_states.rs)

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

**Test:** [`packages/treetime/src/representation/partition/marginal_helpers.rs`](../../packages/treetime/src/representation/partition/marginal_helpers.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/representation/partition/marginal_helpers.rs`](../../packages/treetime/src/representation/partition/marginal_helpers.rs)

| Test                                   | Purpose                                                                    |
| -------------------------------------- | -------------------------------------------------------------------------- |
| `test_propagate_raw_per_site_forward`  | Forward propagation with per-site rates matches individual expQt_with_rate |
| `test_propagate_raw_per_site_backward` | Backward propagation (transpose) with per-site rates matches individual    |

**Test:** [`packages/treetime/src/representation/payload/ancestral.rs`](../../packages/treetime/src/representation/payload/ancestral.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/representation/payload/ancestral.rs`](../../packages/treetime/src/representation/payload/ancestral.rs)

| Test                                                       | Purpose                                              |
| ---------------------------------------------------------- | ---------------------------------------------------- |
| `test_annotate_branch_mutations_formats_1_based_positions` | Mutations formatted with 1-based positions           |
| `test_annotate_branch_mutations_empty_partitions`          | Empty partition list produces no mutations           |
| `test_annotate_branch_mutations_no_mutations_on_edge`      | Edge with no subs produces None                      |
| `test_annotate_branch_mutations_sorts_by_position`         | Multiple mutations sorted by position                |
| `test_annotate_branch_mutations_multi_partition_merge`     | Mutations from multiple partitions merged and sorted |

**Test:** [`packages/treetime/src/representation/payload/discrete.rs`](../../packages/treetime/src/representation/payload/discrete.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/representation/payload/discrete.rs`](../../packages/treetime/src/representation/payload/discrete.rs)

| Test                                         | Purpose                                             |
| -------------------------------------------- | --------------------------------------------------- |
| `test_from_observed_creates_one_hot_profile` | Observed index creates one-hot profile              |
| `test_missing_creates_uniform_profile`       | Missing creates uniform profile                     |
| `test_default_node_data`                     | Default node data has empty profile and zero log_lh |
| `test_default_edge_data`                     | Default edge data has empty arrays and zero log_lh  |

**Test:** [`packages/treetime/src/representation/payload/timetree.rs`](../../packages/treetime/src/representation/payload/timetree.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/representation/payload/timetree.rs`](../../packages/treetime/src/representation/payload/timetree.rs)

| Test                                                           | Purpose                                                         |
| -------------------------------------------------------------- | --------------------------------------------------------------- |
| `test_timetree_annotate_branch_mutations_populates_base_field` | annotate_branch_mutations writes to NodeTimetree.base.mutations |
| `test_timetree_nwk_comments_include_mutations_and_date`        | nwk_comments includes both mutations and date annotations       |
| `test_timetree_nexus_output_includes_mutations_and_date`       | Nexus serialization carries mutations and date annotations      |

### Commands: Prune

**Test:** [`packages/treetime/src/commands/prune/__tests__/test_run.rs`](../../packages/treetime/src/commands/prune/__tests__/test_run.rs)

**Impl:** [`packages/treetime/src/commands/prune/run.rs`](../../packages/treetime/src/commands/prune/run.rs)

| Test                                                                         | Purpose                                                      |
| ---------------------------------------------------------------------------- | ------------------------------------------------------------ |
| `test_prune_nodes_basic`                                                     | Basic zero-length edge pruning                               |
| `test_prune_nodes_with_threshold`                                            | Threshold preserves edges above limit                        |
| `test_prune_nodes_preserves_large_edges`                                     | Large edges unaffected by threshold                          |
| `test_prune_nodes_empty_graph`                                               | Empty graph handled without error                            |
| `test_prune_nodes_handles_none_weights`                                      | None-weight edges handled                                    |
| `test_prune_nodes_preserves_terminal_nodes`                                  | Terminal nodes never collapsed                               |
| `test_prune_nodes_complex_tree`                                              | Multi-level short edge collapse                              |
| `test_prune_nodes_prune_empty_preserves_leaves`                              | Empty-edge leaves preserved                                  |
| `test_prune_nodes_prune_empty_internal_nodes`                                | Empty internal edges collapsed                               |
| `test_prune_nodes_prune_empty_none_mutations`                                | Unknown mutations preserved (not treated as empty)           |
| `test_prune_nodes_prune_empty_simple_leaf_case`                              | Single leaf with no muts preserved                           |
| `test_prune_nodes_combined_prune_short_and_empty`                            | Short and empty pruning combined                             |
| `test_prune_nodes_prune_short_threshold_exact`                               | Exact threshold boundary preserved                           |
| `test_prune_nodes_prune_short_threshold_below`                               | Below-threshold leaves preserved                             |
| `test_prune_nodes_prune_empty_complex_tree`                                  | Complex tree empty-edge collapse                             |
| `test_prune_nodes_prune_both_disabled`                                       | No pruning when both options off                             |
| `test_collapse_sparse_edges_from_leaf_recursive_basic`                       | Recursive leaf path removal                                  |
| `test_collapse_sparse_edges_from_leaf_recursive_stops_at_node_with_children` | Stops at node with remaining children                        |
| `test_collapse_sparse_edges_from_leaf_recursive_stops_at_root`               | Stops at root node                                           |
| `test_collapse_sparse_edges_from_leaf_recursive_invalid_edge_key_errors`     | Invalid edge key errors                                      |
| `test_create_test_edge_num_muts_none_vs_some_zero`                           | None vs Some(0) mutation distinction                         |
| `test_prune_nodes_single_named_leaf`                                         | Single named leaf removal                                    |
| `test_prune_nodes_multiple_named_leaves`                                     | Multiple named leaves removal                                |
| `test_prune_nodes_nonexistent_name_is_noop`                                  | Nonexistent name is no-op                                    |
| `test_prune_nodes_named_internal_node`                                       | Named internal node collapsed                                |
| `test_prune_nodes_mixed_internal_and_leaf_names`                             | Mixed internal and leaf removal                              |
| `test_prune_nodes_empty_names_set_is_noop`                                   | Empty names set is no-op                                     |
| `test_prune_nodes_all_leaves_preserves_root`                                 | Removing all leaves keeps root                               |
| `test_prune_nodes_deep_nested_leaf_removal`                                  | Deep nested leaf removal with collapse                       |
| `test_collapse_edge_mutation_union_non_overlapping`                          | Non-overlapping mutations merged on collapse                 |
| `test_collapse_edge_mutation_union_overlapping_same_mutation`                | Identical mutations deduplicated                             |
| `test_collapse_edge_mutation_union_different_mutations_same_position`        | Different mutations at same position kept                    |
| `test_collapse_edge_mutation_union_multiple_partitions`                      | Per-partition mutation merge                                 |
| `test_collapse_edge_branch_length_sum_both_some`                             | Both branch lengths summed                                   |
| `test_collapse_edge_branch_length_sum_precision`                             | Small branch length sum precision                            |
| `test_collapse_edge_branch_length_none_plus_some`                            | None + Some stays None                                       |
| `test_collapse_edge_branch_length_some_plus_none`                            | Some + None preserves Some                                   |
| `test_collapse_edge_branch_length_both_none`                                 | Both None stays None                                         |
| `test_collapse_edge_compose_non_overlapping`                                 | Non-overlapping subs preserved through edge collapse         |
| `test_collapse_edge_compose_chain`                                           | Chain composition (A->G + G->T = A->T) through edge collapse |
| `test_collapse_edge_compose_cancellation`                                    | Cancellation (A->G + G->A = none) through edge collapse      |
| `test_collapse_edge_compose_multiple_partitions`                             | Composition applied independently per partition              |

**Test:** [`packages/treetime/src/commands/prune/__tests__/test_merge_shared_mutations.rs`](../../packages/treetime/src/commands/prune/__tests__/test_merge_shared_mutations.rs)

**Impl:** [`packages/treetime/src/commands/prune/run.rs`](../../packages/treetime/src/commands/prune/run.rs)

| Test                                                      | Purpose                                               |
| --------------------------------------------------------- | ----------------------------------------------------- |
| `test_merge_no_polytomy`                                  | Binary tree unchanged                                 |
| `test_merge_polytomy_no_shared_mutations`                 | No shared mutations, tree unchanged                   |
| `test_merge_polytomy_two_siblings_share_all_mutations`    | Identical mutation sets merged under new node         |
| `test_merge_polytomy_partial_overlap`                     | Partial overlap: shared moved, unique retained        |
| `test_merge_greedy_picks_best_pair`                       | Pair with most shared mutations selected first        |
| `test_merge_branch_length_adjustment`                     | New edge bl = JC69-corrected distance, child adjusted |
| `test_merge_branch_length_clamp_to_zero`                  | Child bl clamped to zero when smaller than new bl     |
| `test_merge_multiple_partitions`                          | Shared mutations computed per-partition               |
| `test_merge_repeated_until_exhausted`                     | Multiple merge rounds for independent pairs           |
| `test_merge_empty_partitions`                             | No partitions, no merges                              |
| `test_merge_preserves_tree_structure_for_non_polytomies`  | Binary subtrees unaffected                            |
| `test_merge_single_mutation_shared`                       | Minimal case: single shared mutation                  |
| `test_merge_branch_length_jc_correction_differs_from_raw` | Regression guard: JC69 correction applied, not raw p  |
