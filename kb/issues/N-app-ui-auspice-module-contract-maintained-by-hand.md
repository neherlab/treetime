# The Auspice module declarations are checked by hand

The UI embeds Auspice by importing its source modules directly (`auspice/src/...`), because Auspice publishes no library build. `packages/app-ui/tsconfig.json` maps `auspice/src/*` to hand-written declarations in `packages/app-ui/src/auspice/types/`, one `.d.ts` per imported module. Nothing checks the declarations against the Auspice sources, so an Auspice upgrade can change a module's exports or state shape without any compile error.

## Evidence

- The declarations in `packages/app-ui/src/auspice/types/` type the imports in `packages/app-ui/src/auspice/`. `fn auspiceEntries()` in `packages/app-ui/build/auspice-vite.ts` derives the modules that Vite pre-bundles from the same files, so a declared module is always pre-bundled, and TypeScript fails on an import without a declaration
- `packages/app-ui/src/auspice/state.ts` describes the part of the Auspice Redux state the UI reads. Auspice's own types for these modules (`ControlsState`, `TreeState`, `RootState` in `src/reducers/`) never reach the UI because the path mapping replaces them

## Impact

A declaration can diverge from the Auspice module it describes. For example, Auspice types `ScatterVariables` domains as `number[]` and has no `xLabel`/`yLabel` fields, but its scatterplot middleware writes string domains and labels at runtime. The declaration of `ScatterVariables` in `packages/app-ui/src/auspice/state.ts` follows the runtime behavior, and nothing flags a later change on either side.

## Possible changes

- Review the declarations under `packages/app-ui/src/auspice/types/` against the Auspice source at each Auspice upgrade

Using Auspice's own TypeScript types directly would type-check Auspice's `.ts` sources under the UI's strict compiler settings, so the hand-written declarations stay.

## Proposed implementation plan

1. Add a check to the Auspice upgrade procedure (the patch procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) and the developer guide in `docs/en/development/`): for each `.d.ts` under `packages/app-ui/src/auspice/types/`, diff the matching module between the old and the new Auspice tag (`git -C .repos/auspice diff <old>..<new> -- src/<module path>`), and diff `src/reducers/` against the state shape in `packages/app-ui/src/auspice/state.ts`
2. Add a routine TypeScript test in `packages/app-ui/src/auspice/__tests__/` that imports each declared module through Vite's resolution of the real Auspice source and asserts that its runtime export names equal the names that its declaration lists. Oracle: the Auspice source module, which is independent of the hand-written declaration. The list of declared modules comes from the same files that `fn auspiceEntries()` in `packages/app-ui/build/auspice-vite.ts` reads
3. Keep the state-shape fields (such as the string domains and `xLabel`/`yLabel` of `ScatterVariables`) under the review of step 1, because a runtime export check cannot see field types

Validation: `just typecheck` and `just test-ts`; see red by removing one export name from a declaration and checking that the new test fails with that name.
