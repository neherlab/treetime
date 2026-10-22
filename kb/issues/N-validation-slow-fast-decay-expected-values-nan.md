# Validation case `slow_fast_decay` fails on NaN expected values

The case `slow_fast_decay` of the suite `conv-exponential-exponential` in `treetime-validation` fails for every convolution algorithm, with the rectangle sum and with the trapezoid rule, before any convolution result is compared:

```
When finding the peak of the expected values: Undefined ordering between a tested pair of values
```

The error comes from `expected.max()` in `fn compute_peak_metrics()` [packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/peak_metrics.rs#L17](../../packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/peak_metrics.rs#L17): `ndarray-stats` refuses the maximum of an array that holds a NaN. The analytic expected values of the case therefore hold a NaN. The case has no result in the suite, so the role swap it was written to test (rates $a = 0.1$ and $b = 10$ on $[0, 200]$, 200001 points, the swap of `fast_slow_decay`) is not validated.

The case is defined in [packages/treetime-analytical/src/validation/cases/exponential_convolution.rs#L79-L89](../../packages/treetime-analytical/src/validation/cases/exponential_convolution.rs#L79-L89), and the expected values come from `fn exponential_convolution()` [packages/treetime-analytical/src/exponential.rs#L16-L25](../../packages/treetime-analytical/src/exponential.rs#L16-L25), which evaluates $\frac{ab}{a-b}\left(1 - e^{-(a-b)x}\right)e^{-bx}$.

> [!IMPORTANT]
> **Investigation required.** The source of the NaN. Candidate from reading the formula: for $a < b$ the factor $e^{-(a-b)x} = e^{9.9x}$ overflows to infinity above $x \approx 72$, and $e^{-bx} = e^{-10x}$ underflows to 0 above $x \approx 75$, so their product is $\infty \cdot 0 = \text{NaN}$. Evidence to collect: the indices of the NaN values in the expected array of the case, and whether `fast_slow_decay` ($a > b$) avoids them.

## Fix direction

Evaluate the closed form without the overflow, for example as $\frac{ab}{b-a}\left(e^{-ax} - e^{-bx}\right)$ through `exp_m1`, which stays finite for both orders of the rates, and test `fn exponential_convolution()` at large $x$ for $a < b$ and $a > b$.

## Validation

`slow_fast_decay` runs through every metric for every algorithm, and its expected values are finite at every grid point.
