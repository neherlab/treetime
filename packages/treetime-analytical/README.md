# treetime-analytical

Closed-form analytical formulas for probability distributions used in phylogenetic inference. Provides exact solutions for operations that would otherwise require numerical approximation.

## Features

### Gaussian Operations

- `gaussian_pdf` / `gaussian_pdf_grid` - Evaluate normalized Gaussian PDF
- `gaussian_evaluate` - Evaluate Gaussian with custom amplitude on grid
- `gaussian_product` / `gaussian_product_params` - Compute product of N Gaussians (also Gaussian)
- `gaussian_convolution` / `gaussian_convolution_pdf` - Convolve two Gaussians (also Gaussian)

### Exponential Operations

- `exponential_pdf` / `exponential_pdf_grid` - Evaluate exponential PDF
- `exponential_convolution` / `exponential_convolution_grid` - Convolve two exponentials (handles equal/distinct rates)

### Mixed Operations

- `gaussian_exponential_convolution` - Convolve exponential with standard Gaussian (uses complementary error function)

### Validation Test Cases

Pre-defined test cases in `validation::cases` for verifying numerical algorithms against analytical ground truth:

- Gaussian convolution
- Gaussian multiplication
- Gaussian chain multiplication
- Exponential convolution
- Gaussian-exponential convolution

## Usage

```rust
use ndarray::array;
use treetime_analytical::{
    GaussianParams, gaussian_product, gaussian_convolution,
    exponential_convolution_grid,
};

// Product of two Gaussians
let params = vec![
    GaussianParams { mu: 0.0, sigma: 1.0, amplitude: 1.0 },
    GaussianParams { mu: 1.0, sigma: 2.0, amplitude: 1.0 },
];
let grid = array![-5.0, -2.5, 0.0, 2.5, 5.0];
let result = gaussian_product(&params, &grid);
// result.normalized contains the shape, result.log_scale the scale factor

// Convolution of two Gaussians
let a = GaussianParams { mu: 0.0, sigma: 1.0, amplitude: 1.0 };
let b = GaussianParams { mu: 1.0, sigma: 1.0, amplitude: 1.0 };
let conv = gaussian_convolution(&a, &b, &grid);

// Convolution of two exponentials
let exp_conv = exponential_convolution_grid(1.0, 2.0, &grid);
```

## Mathematical Background

### Gaussian Product

The product of N Gaussians is also Gaussian. For Gaussians with means `mu_i` and standard deviations `sigma_i`:

- Combined precision: `1/sigma*^2 = sum(1/sigma_i^2)`
- Combined mean: `mu* = sigma*^2 * sum(mu_i/sigma_i^2)`

### Gaussian Convolution

The convolution of two Gaussians is Gaussian with:

- `mu = mu_1 + mu_2`
- `sigma = sqrt(sigma_1^2 + sigma_2^2)`

### Exponential Convolution

For exponentials with rates `a` and `b`:

- If `a != b`: `(a*b)/(a-b) * (1 - exp(-(a-b)*x)) * exp(-b*x)`
- If `a == b` (limit): `a*b*x * exp(-a*x)`
