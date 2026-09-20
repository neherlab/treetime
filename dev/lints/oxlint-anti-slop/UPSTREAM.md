# Vendored anti-slop Oxlint plugin

TreeTime vendors the `anti-slop` Oxlint plugin instead of depending on an npm package. Upstream ships no official package and documents vendoring as its distribution model.

## Source

- Repository: <https://github.com/dmmulroy/anti-slop>
- Commit: `c44ef22ca116d0ba62a3ff663a0bd13a3f3fa40b`
- Upstream source path: `src/`
- Local destination: `dev/lints/oxlint-anti-slop/` (the upstream `src/` level is flattened into this directory, so the generic entry point is `dev/lints/oxlint-anti-slop/index.ts`)
- License: MIT (`LICENSE`, `Copyright (c) 2026 Dillon Mulroy`), retained verbatim
- Nested vendored code: `vendor/eslint-stylistic/` carries its own MIT `LICENSE` (OpenJS Foundation and ESLint Stylistic contributors) and `UPSTREAM.md`, both retained verbatim

The vendored tree is the complete upstream `src/` tree: the generic plugin (`index.ts`, `rules/`, `shared/`), the Effect plugin (`effect/`), every rule test, and the nested eslint-stylistic vendor. Upstream package-manager files (`package.json`, `pnpm-lock.yaml`, `tsconfig.json`) and the skill-packaging script are not copied; TreeTime type-checks and tests the tree through its own configuration.

## Loading and activation

Oxlint loads the generic plugin from `./dev/lints/oxlint-anti-slop/index.ts`, registered in `oxlint.config.ts` under the plugin name `anti-slop`. The Effect plugin (`anti-slop-effect`, `effect/index.ts`) stays vendored but unregistered because TreeTime has no direct Effect dependency.

TreeTime configuration owns rule severity and path overrides. Every generic rule has one explicit disposition:

### Enabled as errors

- `anti-slop/no-array-filter-map`
- `anti-slop/no-conditional-empty-object-spread`
- `anti-slop/no-known-value-widening`
- `anti-slop/no-module-mocking`
- `anti-slop/no-object-parameters`
- `anti-slop/no-reduce-accumulator-copy`
- `anti-slop/no-reflect-apply`
- `anti-slop/no-reflect-get`
- `anti-slop/no-runtime-typeof` (configured with `allowInTypeGuards: true`)
- `anti-slop/no-shape-in-symbol-names`
- `anti-slop/no-unknown-parameters`
- `anti-slop/no-unknown-returns`
- `anti-slop/no-unknown-type-aliases`
- `anti-slop/no-unsafe-dictionary-type`
- `anti-slop/require-readable-spacing`

`oxc/no-accumulating-spread` is enabled alongside `anti-slop/no-reduce-accumulator-copy`; the two cover growing-accumulator copies from complementary directions.

### Disabled: superseded by a stricter TreeTime rule

- `anti-slop/no-chained-type-assertions`
- `anti-slop/no-widen-then-assert`
- `anti-slop/require-safety-comment-for-type-assertion`

TreeTime keeps `typescript/consistent-type-assertions` with `assertionStyle: "never"`, which rejects every non-const type assertion earlier and more strictly than the three assertion rules above. Enabling them would only add duplicate diagnostics for assertions that are already forbidden.

### Replaced first-party rule

`anti-slop/no-module-mocking` owns the module-mocking contract. The former `treetime/no-module-mocks` rule is removed so one rule owns the policy.

## Local patch

- `rules/require-readable-spacing-cli.test.ts`: the child-process command is changed from `pnpm exec oxlint` to `bun --bun oxlint`. TreeTime runs Oxlint through Bun with the pinned workspace version; `pnpm` is not installed in the container. No other upstream file is modified.

Local functional changes to the vendored tree require an entry in this file and a regression test.

## Style enforcement

The vendored tree is excluded from Oxlint and Oxfmt style enforcement (see the `ignorePatterns` in `oxlint.config.ts` and `oxfmt.config.ts`). It is type-checked with its own `tsconfig.vendor.json`, which mirrors upstream's compiler options rather than TreeTime's stricter first-party flags, so the tree type-checks verbatim. Its rule tests run in the standard custom-rule test path.

## Verification

- Rule tests: `./dev/docker/run just oxlint-test`
- Type check: `./dev/docker/run bun run typecheck:vendor`
- Full read-only gate: `./dev/docker/run just check-all`

## Updating

1. Fetch the target upstream revision into a read-only clone and record its commit.
2. Compare the clone's `src/` against this directory (`diff -r`), classifying every difference as an upstream change or a recorded local patch.
3. Copy the new upstream source, re-apply the local patch above, and update the source commit and any changed rule dispositions in this file.
4. Preserve every license and provenance file verbatim.
5. Run the verification commands and reconcile any new or changed diagnostics.
