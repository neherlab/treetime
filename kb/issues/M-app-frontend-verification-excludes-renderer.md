# Frontend verification excludes the renderer

The desktop `typecheck` script checks only the main-process configuration. The renderer configuration is not invoked and fails directly because CSS ambient declarations are absent. The web test command exits unsuccessfully because no tests exist, and formatting fails in `ResultsPanel.tsx`.

## Evidence

- The desktop `typecheck` script invokes only the default TypeScript project [packages/app-desktop/package.json#L12-L17](../../packages/app-desktop/package.json#L12-L17), while the renderer has a separate configuration [packages/app-desktop/tsconfig.renderer.json#L1-L18](../../packages/app-desktop/tsconfig.renderer.json#L1-L18).
- The web package declares its checks independently [packages/app-web/package.json#L12-L18](../../packages/app-web/package.json#L12-L18), so workspace verification must prove that tests are collected rather than treating an empty suite as coverage.

## Potential solutions

- O1. Make the workspace check invoke both desktop TypeScript projects explicitly.
- O2. Use a composite TypeScript project with references to main and renderer configs.

## Recommendation

Make declared workspace checks cover the renderer, add the required Vite/CSS ambient declarations, add tests for transport and request construction, and keep the formatter clean.

## Related issues

- [H-app-transport-contracts-diverge-across-clients.md](H-app-transport-contracts-diverge-across-clients.md)
- [H-app-ui-displays-synthetic-results-and-ignores-inputs.md](H-app-ui-displays-synthetic-results-and-ignores-inputs.md)
