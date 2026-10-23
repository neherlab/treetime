# The Auspice module declarations are checked by hand

The UI embeds Auspice by importing its source modules directly (`auspice/src/...`), because Auspice publishes no library build. `packages/app-ui/tsconfig.json` maps `auspice/src/*` to hand-written declarations in `packages/app-ui/src/auspice/types/`, one `.d.ts` per imported module. Nothing checks the declarations against the Auspice sources, so an Auspice upgrade can change a module's exports or state shape without any compile error.

## Evidence

- The declarations in `packages/app-ui/src/auspice/types/` type the imports in `packages/app-ui/src/auspice/`. `fn auspiceEntries()` in `packages/app-ui/build/auspice-vite.ts` derives the modules that Vite pre-bundles from the same files, so a declared module is always pre-bundled, and TypeScript fails on an import without a declaration
- `packages/app-ui/src/auspice/state.ts` describes the part of the Auspice Redux state the UI reads. Auspice's own types for these modules (`ControlsState`, `TreeState`, `RootState` in `src/reducers/`) never reach the UI because the path mapping replaces them

## Impact

A declaration can diverge from the Auspice module it describes. For example, Auspice types `ScatterVariables` domains as `number[]` and has no `xLabel`/`yLabel` fields, but its scatterplot middleware writes string domains and labels at runtime. The declaration here follows the runtime behavior, and nothing flags a later change on either side.

## Possible changes

- Review the declarations under `packages/app-ui/src/auspice/types/` against the Auspice source at each Auspice upgrade

Using Auspice's own TypeScript types directly would type-check Auspice's `.ts` sources under the UI's strict compiler settings, so the hand-written declarations stay.
