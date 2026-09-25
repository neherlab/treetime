# The effect lint rules report effects of earlier files at the end of later files

## Summary

`web/no-state-in-effect-initializer` collects effects and state setters in a state object that `makeEffectState()` creates once per `createOnce` call, which is once per linter process. Nothing resets it between files, so after one file contains a violating mount effect, the rule reports that effect again at `Program:exit` of every later file, pointing at unrelated files and line numbers.

## Evidence

- `fn makeEffectState()` builds `{ setterNames, effects }` once and the visitors only append to it [`dev/lints/oxlint/rules/effects.ts#L38`](../../dev/lints/oxlint/rules/effects.ts#L38)
- `noStateInEffectInitializerRule` iterates `state.effects` at `Program:exit` without clearing them [`dev/lints/oxlint/rules/no-state-in-effect-initializer.ts`](../../dev/lints/oxlint/rules/no-state-in-effect-initializer.ts)
- A single mount effect that set state in a `ResizeObserver` callback produced 38 reports, most of them at the last line of files without any effect, such as `packages/app-ui/src/ui/fonts.ts` and `packages/app-web/src/env.d.ts`

## Impact

- One real finding floods the lint output with false locations, which hides the file that needs the fix
- Any other rule that uses `makeEffectState()` shares the defect
