# treetime-analytical

Closed-form analytical formulas for probability distributions used in phylogenetic inference. Provides exact solutions for operations that would otherwise require numerical approximation.

## Operations

### Gaussian

Core parameter type: `GaussianParams { mu, sigma, amplitude }`.

- `gaussian_pdf` / `gaussian_pdf_grid` - Normalized Gaussian PDF: `exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))`
- `gaussian_evaluate` - Evaluate Gaussian with custom amplitude on grid
- `gaussian_product` / `gaussian_product_params` - Product of N Gaussians (also Gaussian), returns `ScaledArray` with normalized shape and log-scale factor
- `gaussian_convolution` - Convolve two Gaussians with arbitrary amplitude
- `gaussian_convolution_pdf` / `gaussian_convolution_pdf_grid` - Convolve two normalized Gaussian PDFs

```rust
use ndarray::array;
use treetime_analytical::{GaussianParams, gaussian_product, gaussian_convolution};

let params = vec![
    GaussianParams { mu: 0.0, sigma: 1.0, amplitude: 1.0 },
    GaussianParams { mu: 1.0, sigma: 2.0, amplitude: 1.0 },
];
let grid = array![-5.0, -2.5, 0.0, 2.5, 5.0];

// Product of two Gaussians
let result = gaussian_product(&params, &grid);
// result.normalized contains the shape, result.log_scale the scale factor

// Convolution of two Gaussians
let a = GaussianParams { mu: 0.0, sigma: 1.0, amplitude: 1.0 };
let b = GaussianParams { mu: 1.0, sigma: 1.0, amplitude: 1.0 };
let conv = gaussian_convolution(&a, &b, &grid);
```

### Exponential

- `exponential_pdf` / `exponential_pdf_grid` - Exponential PDF: `rate * exp(-rate * x)` for `x >= 0`
- `exponential_convolution` / `exponential_convolution_grid` - Convolve two exponentials, handling both distinct and equal rate cases

```rust
use ndarray::array;
use treetime_analytical::{exponential_pdf, exponential_convolution_grid};

let grid = array![0.0, 0.5, 1.0, 1.5, 2.0];

let pdf = exponential_pdf(2.0, 1.0);
let conv = exponential_convolution_grid(1.0, 2.0, &grid);
```

### Gaussian-exponential convolution

- `gaussian_exponential_convolution` / `gaussian_exponential_convolution_grid` - Convolve exponential with standard Gaussian using the complementary error function

```rust
use ndarray::array;
use treetime_analytical::gaussian_exponential_convolution_grid;

let grid = array![-2.0, 0.0, 2.0, 5.0];
let result = gaussian_exponential_convolution_grid(0.5, &grid);
```

## Validation test cases

Pre-defined test cases in `validation::cases` for verifying numerical algorithms against analytical ground truth. Each case specifies distribution parameters, grid domain, and stress characteristics.

- `GAUSSIAN_CONVOLUTION_CASES` - Gaussian convolution with varying sigma ratios and offsets
- `GAUSSIAN_PAIRWISE_MULTIPLICATION_CASES` - Pairwise Gaussian multiplication with scale extremes
- `get_gaussian_chain_multiplication_cases()` - Chained Gaussian multiplication sequences
- `EXPONENTIAL_CONVOLUTION_CASES` - Exponential convolution with distinct and near-equal rates
- `GAUSSIAN_EXPONENTIAL_CASES` - Mixed Gaussian-exponential convolution

## Mathematical background

### Gaussian product

The product of N Gaussians is Gaussian with combined precision `1/sigma*^2 = sum(1/sigma_i^2)` and mean `mu* = sigma*^2 * sum(mu_i / sigma_i^2)`.

### Gaussian convolution

The convolution of two Gaussians has `mu = mu_1 + mu_2` and `sigma = sqrt(sigma_1^2 + sigma_2^2)`.

### Exponential convolution

For exponentials with rates `a` and `b`:

- `a != b`: `(a*b)/(a-b) * (1 - exp(-(a-b)*x)) * exp(-b*x)`
- `a == b` (limit): `a*b*x * exp(-a*x)`

### Gaussian-exponential convolution

For exponential with rate `a` and standard Gaussian: `0.5 * a * exp(-x*a + 0.5*a^2) * erfc((a - x) / sqrt(2))`.
