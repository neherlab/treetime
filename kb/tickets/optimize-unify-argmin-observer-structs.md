# Unify argmin observer structs into single parameterized type

Three optimization sites define near-identical `Observe<I>` implementations differing only in label and early-iteration threshold.

## Current state

`TcOptimizationObserver` (coalescent), `BrentObserver` (clock/find_best_root), and `GoldenSectionObserver` (clock/find_best_root) are separate unit structs with identical trait bounds and body structure.

## Target state

A single `OptimizationObserver { label: &'static str, early_threshold: u64 }` struct with one `Observe` implementation. All three sites instantiate with their specific label and threshold.

## Implementation

1. Define `OptimizationObserver` in a shared location (e.g., `optimize/observer.rs` or `optimize/mod.rs`)
2. Implement `Observe<I>` once with parameterized label and threshold
3. Replace `TcOptimizationObserver` in `coalescent/optimize_tc.rs:139-159` with `OptimizationObserver { label: "Tc optimization", early_threshold: 3 }`
4. Replace `BrentObserver` in `clock/find_best_root/method_brent.rs:13-33` with `OptimizationObserver { label: "Brent", early_threshold: 5 }`
5. Replace `GoldenSectionObserver` in `clock/find_best_root/method_golden_section.rs:13-33` with `OptimizationObserver { label: "Golden Section", early_threshold: 5 }`
6. Delete the three original structs

## Related issues

Source: `kb/issues/N-optimize-argmin-observer-boilerplate.md`
