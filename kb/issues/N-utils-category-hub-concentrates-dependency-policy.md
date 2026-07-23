# treetime-utils category hub concentrates dependency policy

`treetime-utils` groups arrays, collections, dates, errors, formatting, initialization, intervals, I/O, iterators, synchronization, and testing under one crate root [`packages/treetime-utils/src/lib.rs`](../../packages/treetime-utils/src/lib.rs). Its leaf modules are cohesive, but unrelated workspace crates share the same stable crate-level dependency.

## Risk

The crate has broad afferent coupling: changes to root dependencies, features, initialization, or layering policy can affect consumers that use unrelated leaf utilities. Domain-specific behavior added for convenience would make the most widely shared workspace node depend on one domain and invert ownership.

This does not justify splitting a stable crate solely by module count. The actionable defect is the absence of a boundary rule preventing new cross-domain coupling in the shared hub.

## Required policy

- Preserve leaf-module APIs and keep the root limited to module declarations.
- Place domain-specific utilities in their owning crate.
- Do not add cross-leaf mutable state or initialization dependencies.
- If unrelated consumers begin to co-change because of features or dependencies, split along existing top-level module boundaries with semantic impact analysis.

## Validation

- Dependency review accompanies every new `treetime-utils` dependency or feature.
- Workspace dependency analysis identifies which leaf module each consumer imports.
- Domain crates do not appear in the utility crate's dependency graph.
- A split is proposed only with observed co-change or dependency pressure, not module count alone.
