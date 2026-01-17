# treetime-ops

Numerical operations for probability distributions on uniform grids.

## Overview

Provide convolution and multiplication algorithms for discrete probability distributions represented as `ndarray::Array1<f64>`. These operations are fundamental to phylogenetic likelihood calculations in Treetime.

## Operations

### Convolution

Convolve two distributions on a uniform grid with spacing `dx`. Three implementations available:

- `RiemannConvolve` - Simple O(n\*m) Riemann sum, suitable for small arrays
- `NdarrayConvolve` - Direct convolution via ndarray-conv
- `FftConvolve` - FFT-based O(n log n) convolution, best for large arrays

```rust
use treetime_ops::{convolve_fft, FftConvolve, ConvolveAlgo};
use ndarray::array;

let dx = 0.1;
let f = array![1.0, 2.0, 3.0, 2.0, 1.0];
let g = array![1.0, 1.0, 1.0];

// Function form
let result = convolve_fft(dx, &f, &g).unwrap();

// Trait form
let algo = FftConvolve;
let result = algo.convolve(dx, &f, &g).unwrap();
```

### Multiplication

Multiply multiple distributions element-wise with log-scale tracking to prevent underflow. Returns `ScaledArray` containing normalized values and accumulated log-scale.

- `PointwiseMultiply` - Naive multiplication, can underflow with many factors
- `LogScaleMultiply` - Lazy normalization, normalizes only when approaching underflow
- `AggressiveMultiply` - Normalizes after every pairwise multiplication

```rust
use treetime_ops::{multiply_many_lazy_normalize, LogScaleMultiply, MultiplyAlgo};
use ndarray::array;

let a = array![0.001, 0.002, 0.001];
let b = array![0.5, 1.0, 0.5];
let c = array![0.1, 0.2, 0.1];

// Function form
let result = multiply_many_lazy_normalize(&[&a, &b, &c]);
// result.normalized contains values with max = 1.0
// result.log_scale contains ln(scale_factor)
// Full result = result.normalized * exp(result.log_scale)

// Trait form
let algo = LogScaleMultiply;
let result = algo.multiply_many(&[&a, &b, &c]);
```

## Types

### ScaledArray

Container for normalized distribution and its log-scale factor:

```rust
pub struct ScaledArray {
    pub normalized: Array1<f64>,  // Values normalized so max = 1.0
    pub log_scale: f64,           // ln(scale_factor)
}
```

Reconstruct full values as `normalized * exp(log_scale)`.

### Traits

- `ConvolveAlgo` - Trait for convolution algorithm implementations
- `MultiplyAlgo` - Trait for multiplication algorithm implementations

Both traits are `Send + Sync` for use in parallel computations.
