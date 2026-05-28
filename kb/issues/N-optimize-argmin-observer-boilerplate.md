# Argmin observer structs duplicated across three optimization sites

Three near-identical `Observe<I>` implementations exist for logging argmin optimization progress. Each defines a unit struct with identical trait bounds, identical body (log on multiples of 10 and early iterations), differing only in label string and early-iteration threshold (3 vs 5).

## Locations

- `packages/treetime/src/coalescent/optimize_tc.rs:139-159:` `TcOptimizationObserver` - logs if `iter % 10 == 0 || iter <= 3`, label `"Tc optimization"`
- `packages/treetime/src/clock/find_best_root/method_brent.rs:13-33:` `BrentObserver` - logs if `iter % 10 == 0 || iter <= 5`, label `"Brent"`
- `packages/treetime/src/clock/find_best_root/method_golden_section.rs:13-33:` `GoldenSectionObserver` - logs if `iter % 10 == 0 || iter <= 5`, label `"Golden Section"`

All share identical trait bounds:

```rust
impl<I> Observe<I> for XxxObserver
where
  I: State,
  <I as State>::Param: std::fmt::Debug,
  <I as State>::Float: std::fmt::LowerExp,
```

## Impact

Pure maintainability. Adding a new optimizer requires copying the boilerplate.

## Action

Extract a single `OptimizationObserver { label: &'static str, early_threshold: u64 }` struct with one `Observe` implementation. All three sites instantiate with their label and threshold.
