# Reusable library selects the native linear-algebra backend

> [!IMPORTANT]
> **Decision required.** The supported application, target, and backend matrix is unapproved, and so is the backend selection mechanism: one backend in each final application crate (O1), or documented top-level backend features (O2). See "Potential solutions".

The reusable library dependency graph enables `openblas-system`. Cargo features are additive, so downstream applications cannot disable it when selecting a different `ndarray-linalg` backend. Multiple native backend features can then be active together.

## Evidence

Workspace and crate manifests select `ndarray-linalg` features in reusable dependency declarations [Cargo.toml](../../Cargo.toml) [packages/treetime/Cargo.toml](../../packages/treetime/Cargo.toml). Cargo feature unification makes the resulting backend selection a property of the complete application graph.

## Potential solutions

- O1. Select exactly one backend in final application crates.
- O2. Expose mutually documented top-level backend features and require applications to choose one.

## Recommendation

Keep reusable crates backend-neutral. Select exactly one backend in each final application crate or through mutually documented top-level features. Validate feature graphs for every supported build target.

Implementation waits until the supported application/target/backend matrix is approved; changing backend features without that matrix can make supported targets unbuildable or select an unintended native library.

## Related issues
