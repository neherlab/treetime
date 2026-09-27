# The Auspice modules used by the UI are listed and typed by hand in several places

The UI embeds Auspice by importing its source modules directly (`auspice/src/...`), because Auspice publishes no library build. The set of imported modules and their types are kept by hand in three places that nothing checks against each other. An Auspice upgrade can change a module's exports or state shape without any compile error.

## Evidence

- `const AUSPICE_ENTRIES` in `packages/app-ui/build/auspice-vite.ts` lists every imported Auspice module so that Vite pre-bundles it (`optimizeDeps.include`)
- `packages/app-ui/tsconfig.json` maps `auspice/src/*` to hand-written declarations in `packages/app-ui/src/auspice/types/`, one `.d.ts` per imported module
- The imports themselves are in `packages/app-ui/src/auspice/` (`AuspiceTree.tsx`, `store.ts`, `store-hooks.ts`, `i18n.ts`)
- `packages/app-ui/src/auspice/state.ts` describes the part of the Auspice Redux state the UI reads. Auspice's own types for these modules (`ControlsState`, `TreeState`, `RootState` in `src/reducers/`) never reach the UI because the path mapping replaces them

## Impact

- A module imported without an `AUSPICE_ENTRIES` entry still works, but Vite finds it only while serving the page, then re-bundles dependencies and reloads the page
- A declaration can diverge from the Auspice module it describes. For example, Auspice 3.0.0 types `ScatterVariables` domains as `number[]` and has no `xLabel`/`yLabel` fields, but its scatterplot middleware writes string domains and labels at runtime. The declaration here follows the runtime behavior, and nothing flags a later change on either side

## Possible changes

- Derive `AUSPICE_ENTRIES` from the declaration files, or add a test that compares the entry list, the declaration files, and the `auspice/src/` imports
- Review the declarations under `packages/app-ui/src/auspice/types/` against the Auspice source at each Auspice upgrade

Using Auspice's own TypeScript types directly would type-check Auspice's `.ts` sources under the UI's strict compiler settings, so the hand-written declarations stay.
