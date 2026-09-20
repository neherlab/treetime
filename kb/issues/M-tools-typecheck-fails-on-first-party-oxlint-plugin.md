# `typecheck:tools` fails on the first-party Oxlint plugin

## Symptom and reproduction

`bun run typecheck:tools` (the `typescript-config` step of `just check` and `just check-all`) reports type errors and exits non-zero.

```bash
./dev/docker/run bash -c "bun run typecheck:tools"
```

The dominant error is `TS5097: An import path can only end with a '.ts' extension when 'allowImportingTsExtensions' is enabled`, emitted for the many `./foo.ts` relative imports across `dev/lints/oxlint/`. `tsconfig.tools.json` does not set `allowImportingTsExtensions`, so every `.ts` import in the first-party plugin is rejected.

## Impact and scope

- `typecheck:tools` covers the first-party Oxlint plugin (`dev/lints/oxlint/**`) and the root config files (`oxlint.config.ts`, `oxfmt.config.ts`, `vitest.config.ts`).
- The gate has never passed: the plugin source uses `.ts` import specifiers (required by the Node 24 rule-test runner and by Oxlint's `jsPlugins` loader), which `tsc` rejects without `allowImportingTsExtensions`.
- Runtime is unaffected: the rules load and their tests pass under Node 24 type-stripping. This is a type-checking-config defect, not a runtime defect.

## Root cause

Two independent causes:

1. `tsconfig.tools.json` lacks `allowImportingTsExtensions: true`. Adding it (the config has `noEmit: true`, so it is permitted) removes every `TS5097`.
2. After that, latent strict-mode errors remain in the first-party plugin and in `oxlint.config.ts`:
   - `strictNullChecks` / `noUncheckedIndexedAccess` gaps (`possibly undefined` on indexed AST reads in `no-uppercase-test-title.ts`, `no-exec-shell-string.ts`, `no-event-handler-effect.ts`, `no-fake-success.ts`).
   - `@oxlint/plugins` 1.82 ESTree type drift: `ESTree.Identifier` and `ESTree.FunctionExpression` are no longer namespaced members (`ast.ts`, `effects.ts`, `no-vague-identifiers.ts`); the identifier node is `ESTree.IdentifierReference`/`ESTree.IdentifierName` and the function-expression node must be reached with `Extract<ESTree.Node, { type: "FunctionExpression" }>`.
   - `exactOptionalPropertyTypes` incompatibility between the `no-restricted-imports` override literals and Oxlint's `OxlintOverride`/`DummyRuleMap` types in `oxlint.config.ts`.

## Fix approach

1. Add `allowImportingTsExtensions: true` to `tsconfig.tools.json`.
2. Add the missing null guards.
3. Replace the drifted ESTree type references with the current 1.82 names (introduce a shared `Extract<ESTree.Node, { type: "Identifier" }>` alias for identifier nodes).
4. Type the `no-restricted-imports` overrides so they satisfy `OxlintOverride` under `exactOptionalPropertyTypes`.

Then `typecheck:tools` should pass and can be added to continuous integration alongside the vendored-plugin type check.

## Note

The vendored anti-slop plugin (`dev/lints/oxlint-anti-slop/`) type-checks cleanly under its own `tsconfig.vendor.json` and is unaffected by this issue.
