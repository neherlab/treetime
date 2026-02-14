# treetime-utils

General-purpose utilities shared across the Treetime project. Provides error handling macros, ndarray extensions, compression, datetime parsing, formatting, collections, and testing support.

## Modules

### Error Handling (`error`)

Error macros built on `eyre`:

- `make_error!()` - return `Result::Err` with formatted message
- `make_report!()` - create error report for `.map_err()`
- `make_internal_error!()` / `make_internal_report!()` - internal bug errors with developer contact info
- `report_to_string()` - flatten error chain into a single string

### Array (`array`)

#### Ndarray Extensions (`array::ndarray`)

Array operations extending `ndarray`:

- `outer()` - outer product of two vectors
- `argmin_axis()` / `argmax_axis()` - index of min/max along axis
- `min_axis()` / `max_axis()` - reduce along axis
- `minimum()` / `maximum()` / `maximum_scalar()` - element-wise min/max
- `clamp()` / `clamp_min()` / `clamp_max()` - element clamping
- `cumsum_axis()` / `product_axis()` - cumulative sum, product reduction
- `exp()` / `log()` - element-wise math
- `choose1()` / `choose2()` - index-based selection
- `stack_owned()` - stack owned arrays (unlike `ndarray::stack` which requires views)
- `reverse()` / `reverse_inplace()` - array reversal
- `has_uniform_spacing()` - check grid uniformity (ULP-based comparison)
- `to_col()` / `to_row()` - reshape 1D array to column/row matrix
- `ndarray_pad_zeros_right()` - pad or truncate array to target length
- `ndarray_uniform_grid()` - create uniformly spaced grid
- `random()` - generate random array with uniform distribution

#### Serde Helpers (`array::serde`)

Custom serialization for `ndarray` and related types:

- `array1_as_vec` / `array1_from_vec` - serialize `Array1` as JSON array
- `indexmap_array1_from_map` - deserialize `IndexMap<String, Array1<f64>>` from JSON object
- `skip_serializing_if_false` - omit boolean fields when false

### Collections (`collections`)

- `container` - `count_occurrences()`, `get_exactly_one()`, `get_one()`, `minmax()`
- `vec` - `vec_of_owned![]`, `vec_u8![]` macros
- `manyzip` - `Manyzip` iterator for zipping a dynamic number of iterators

### DateTime (`datetime`)

Date and time parsing for phylogenetic metadata:

- `year_fraction` - bidirectional conversion between `DateTime<Utc>` and decimal year fractions (e.g. `2024.679`)
- `date_range` - date range type with begin/end bounds
- `parse_date` / `parse_datetime` - date string parsing
- `parse_uncertain_date` - dates with uncertainty (ranges, partial dates like `2024-XX-XX`)
- `format_to_regex` - convert date format strings to regex patterns
- `options` - date parsing configuration
- `nil` - sentinel/nil date handling

### Formatting (`fmt`)

- `float` - `FloatFormatExt` trait and functions for formatting floats to significant/decimal digits with trailing zero trimming
- `string` - `o!()` ownership macro, `quote()`/`quote_single()`, `truncate()` with directional truncation (left/right/middle) and optional ellipsis

### Initialization (`init`)

- `global` - `global_init()` for `color_eyre` error hooks and `setup_logger()` for colored log output
- `clap_styles` - CLI styling for `clap` argument parser
- `openblas` - OpenBLAS runtime configuration query (`get_openblas_info()`)

### I/O (`io`)

- `compression` - transparent compression/decompression (gzip, bzip2, xz, zstd) with `Compressor`/`Decompressor` streaming wrappers and format detection from file extension
- `file` - buffered file reading/writing with automatic decompression, stdin/stdout support (`open_file_or_stdin()`, `create_file_or_stdout()`)
- `fs` - filesystem helpers: `ensure_dir()`, `absolute_path()`, `extension()`, `has_extension()`, `add_extension()`, `read_file_to_string()`
- `console` - TTY detection (`is_tty()`)

### Iterator (`iterator`)

- `mean_by_key` - `MeanByKey` trait adding `.mean_by_key()` to all iterators
- `union` - `iterator_union()` for sorted, deduplicated merge of two iterators

### Sync (`sync`)

- `mutex` - `unwrap_arc_rwlock()` for extracting inner value from `Arc<RwLock<T>>`
- `random` - seeded RNG creation with `Isaac64Rng`, `random_choice()`, `random_remove()`, `random_pop()`, `random_sequence()`

### Testing (`testing`)

Test assertion macros wrapping `pretty_assertions` and `approx`:

- `pretty_assert_eq!()` - formatted equality assertion with newline normalization
- `pretty_assert_ulps_eq!()` - ULP-based float comparison with pretty diff output
- `pretty_assert_abs_diff_eq!()` - absolute difference float comparison
- `assert_error!()` - assert error message content

## Memory Allocator

On Linux (x86_64/aarch64, glibc/musl), jemalloc is used as the global allocator for improved performance.

## Features

- `png` - enable PNG output for plots (adds `image` and `plotters` bitmap backend)
