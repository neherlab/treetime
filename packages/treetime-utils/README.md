# treetime-utils

Shared utilities for the Treetime project. Provide common functionality for error handling, ndarray operations, compression, datetime parsing, random number generation, and testing support.

## Modules

### Error Handling (`error`)

Error macros built on `eyre`:

- `make_error!()` - return `Result::Err` with formatted message
- `make_report!()` - create error report for `.map_err()`
- `make_internal_error!()` / `make_internal_report!()` - internal bug errors with developer contact info

### Ndarray Extensions (`ndarray`)

Array operations extending `ndarray`:

- `outer()` - outer product of two vectors
- `argmin_axis()` / `argmax_axis()` - index of min/max along axis
- `min_axis()` / `max_axis()` - reduce along axis
- `minimum()` / `maximum()` - element-wise min/max
- `clamp()` / `clamp_min()` / `clamp_max()` - element clamping
- `cumsum_axis()` / `product_axis()` - cumulative sum, product reduction
- `exp()` / `log()` - element-wise math
- `choose1()` / `choose2()` - index-based selection
- `stack_owned()` - stack owned arrays
- `reverse()` / `reverse_inplace()` - array reversal
- `has_uniform_spacing()` - check grid uniformity

### Compression (`compression`)

Transparent compression/decompression based on file extension:

- Supported formats: gzip, bzip2, xz, zstd
- `Compressor` / `Decompressor` - streaming compression wrappers
- `guess_compression_from_filepath()` - detect format from extension

### DateTime (`datetime`)

Date and time parsing for phylogenetic metadata:

- Parse dates with uncertainty (ranges, partial dates)
- Year fraction conversion for molecular clock analysis
- Date range handling

### Interval (`interval`)

Range operations for numeric intervals:

- Union, intersection, difference, complement
- Range properties and queries

### Random (`random`)

Random number generation utilities:

- Seeded RNG creation with `Isaac64Rng`
- `random_choice()` / `random_remove()` / `random_pop()` - selection from collections

### Testing (`assert`)

Test assertion macros:

- `pretty_assert_eq!()` - formatted equality assertion
- `pretty_assert_ulps_eq!()` - ULP-based float comparison
- `pretty_assert_abs_diff_eq!()` - absolute difference comparison
- `assert_error!()` - error message assertion

### Serde Helpers (`serde`)

Custom serialization for `ndarray` types:

- `array1_as_vec` / `array1_from_vec` - serialize `Array1` as JSON array

### Other

- `global_init` - logger and error hook setup
- `clap_styles` - CLI styling
- `console` - terminal utilities
- `fs` / `file` - filesystem helpers
- `string` / `vec` / `container` - collection utilities
- `iterator` - iterator extensions (e.g., `mean_by_key`)
- `manyzip` - multi-iterator zip
- `mutex` - parking_lot mutex helpers
- `openblas` - OpenBLAS configuration

## Memory Allocator

On Linux (x86_64/aarch64, glibc/musl), jemalloc is used as the global allocator for improved performance.

## Features

- `png` - enable PNG output for plots (adds `image` and `plotters` bitmap backend)
