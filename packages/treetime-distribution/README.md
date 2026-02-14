# treetime-distribution

Probability distributions for phylogenetic time inference in Treetime.

## Overview

Represent and manipulate 1D probability distributions over time. Distributions arise during molecular-clock inference, where each node carries a probability distribution over its divergence time. The crate supports arithmetic operations (multiplication, division, convolution, subtraction, negation) needed to propagate these distributions along the tree.

## Distribution Variants

`Distribution<Y>` is an enum generic over a y-axis policy (`Plain` or `NegLog`):

| Variant    | Type                   | Description                                         |
| ---------- | ---------------------- | --------------------------------------------------- |
| `Empty`    | -                      | No information                                      |
| `Point`    | `DistributionPoint`    | Dirac delta at a single time                        |
| `Range`    | `DistributionRange`    | Constant amplitude over an interval                 |
| `Function` | `DistributionFunction` | Tabulated values on a uniform grid (interpolatable) |
| `Formula`  | `DistributionFormula`  | Lazy evaluation via a closure                       |

Type aliases `DistributionPlain` and `DistributionNegLog` fix the policy parameter.

## Y-Axis Policies

The `YAxisPolicy` trait controls how probability values are stored and combined:

- `Plain` - Direct probability values. Multiplication is `a * b`, identity is `1.0`.
- `NegLog` - Negative log probabilities (`-ln(p)`). Multiplication becomes addition, identity is `0.0`.

Marker traits gate operations that are valid only in certain representations:

- `SupportsConvolution` - Implemented for `Plain` only (convolution requires integral sums)
- `SupportsSubtraction` - Implemented for `Plain` only

Conversions between representations: `Distribution<Plain>::to_neglog()` and `Distribution<NegLog>::to_plain()`.

## Operations

All operations are free functions that accept `&Distribution<Y>` arguments:

| Function                                | Description                                     |
| --------------------------------------- | ----------------------------------------------- |
| `distribution_multiplication`           | Pointwise product of two distributions          |
| `distribution_division`                 | Pointwise quotient                              |
| `distribution_convolution`              | Convolution (Plain only)                        |
| `distribution_subtraction`              | Pointwise difference (Plain only)               |
| `distribution_negation`                 | Reflect across time axis: f(x) -> f(-x)         |
| `distribution_scalar_multiplication`    | Scale amplitude by a scalar                     |
| `distribution_map`                      | Apply an arbitrary function to amplitude values |
| `distribution_time_bounds_union`        | Smallest interval covering both distributions   |
| `distribution_time_bounds_intersection` | Overlap interval of two distributions           |

## ScaledDistribution

Decompose a distribution into a normalized shape and a log-scale factor to prevent numerical underflow during repeated multiplications:

```
P(x) = exp(log_scale) * inner(x)     where max(inner) = 1.0
```

Scaled variants of multiplication, division, and convolution operate on `ScaledDistribution` and propagate the log-scale factor through the computation.

## Dependencies

- `treetime-grid` - Uniform grid and interpolation
- `treetime-ops` - Low-level convolution and multiplication algorithms
- `treetime-utils` - Error macros
- `ndarray` - Array storage and arithmetic
