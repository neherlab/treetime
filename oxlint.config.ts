import { readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";

import { defineConfig } from "oxlint";
import type { OxlintOverride } from "oxlint";
import * as z from "zod";

const PACKAGES_DIR = join(import.meta.dirname, "packages");

const WORKSPACE_SCOPE = "@neherlab/";

const IMPORT_BOUNDARY_PATTERNS = [
  {
    group: ["@neherlab/*/src/*", "@neherlab/*/src/**", "@neherlab/*/dist/*", "@neherlab/*/dist/**"],
    message: "Import a workspace package through its entry point, not a deep internal path.",
  },
  {
    group: ["../../*", "../../**"],
    message: "A relative path must not reach into a sibling package. Import it by its @neherlab/* name.",
  },
];

const PACKAGE_GRAPH_MESSAGE =
  "This package may import only the workspace packages it declares. Add the dependency to its package.json, or route through an allowed package.";

interface WorkspacePackage {
  dir: string;
  name: string;
  dependencies: string[];
}

const packageManifestSchema = z.object({
  name: z.string().optional(),
  dependencies: z.record(z.string(), z.string()).optional(),
  devDependencies: z.record(z.string(), z.string()).optional(),
  peerDependencies: z.record(z.string(), z.string()).optional(),
});

function packageBoundaryOverrides(): OxlintOverride[] {
  const packages = workspacePackages();
  const names = packages.map((entry) => entry.name);

  return packages.flatMap((entry) => {
    const allowed = new Set([entry.name, ...entry.dependencies]);
    const banned = names.filter((other) => !allowed.has(other));

    if (banned.length === 0) {
      return [];
    }

    return [
      {
        files: [`packages/${entry.dir}/**`],
        rules: {
          "no-restricted-imports": [
            "error",
            { patterns: [...IMPORT_BOUNDARY_PATTERNS, { group: banned, message: PACKAGE_GRAPH_MESSAGE }] },
          ],
        },
      },
    ];
  });
}

function workspacePackages(): WorkspacePackage[] {
  const packages: WorkspacePackage[] = [];

  for (const dir of readdirSync(PACKAGES_DIR)) {
    let raw: string;

    try {
      raw = readFileSync(join(PACKAGES_DIR, dir, "package.json"), "utf8");
    } catch {
      continue;
    }

    const manifest = packageManifestSchema.parse(JSON.parse(raw));
    const name = manifest.name;

    if (name === undefined || !name.startsWith(WORKSPACE_SCOPE)) {
      continue;
    }

    const declared = {
      ...manifest.dependencies,
      ...manifest.devDependencies,
      ...manifest.peerDependencies,
    };

    const dependencies = Object.keys(declared).filter((key) => key.startsWith(WORKSPACE_SCOPE));
    packages.push({ dir, name, dependencies });
  }

  return packages;
}

export default defineConfig({
  ignorePatterns: [
    "dist",
    "dist-electron",
    "node_modules",
    "bun.lock",
    "*.tsbuildinfo",
    ".turbo",
    ".build",
    "coverage",
    "packages/app-contracts/src/generated",
    "dev/lints/oxlint-anti-slop",
  ],

  plugins: [
    "eslint",
    "typescript",
    "oxc",
    "unicorn",
    "import",
    "promise",
    "node",
    "jsx-a11y",
    "react",
    "react-perf",
    "vitest",
  ],

  jsPlugins: [
    "eslint-plugin-no-comments",
    "eslint-plugin-sonarjs",
    "./dev/lints/oxlint/index.ts",
    "./dev/lints/oxlint/web.ts",
    "./dev/lints/oxlint-anti-slop/index.ts",
  ],

  categories: {
    correctness: "warn",
    suspicious: "warn",
    pedantic: "off",
    perf: "warn",
    style: "off",
  },

  options: {
    typeAware: true,
    denyWarnings: true,
    reportUnusedDisableDirectives: "error",
    respectEslintDisableDirectives: false,
  },

  rules: {
    "no-comments/disallowComments": [
      "error",
      { allow: ["oxlint-", "@ts-expect-error", "\\* @internal", "/ <reference", "/usr/bin/env"] },
    ],

    "sonarjs/class-name": "error",
    "sonarjs/function-name": ["error", { format: "^[_a-z][a-zA-Z0-9]*$|^[A-Z][a-zA-Z0-9]*$" }],
    "sonarjs/variable-name": "error",
    "sonarjs/no-all-duplicated-branches": "error",
    "sonarjs/no-collapsible-if": "error",
    "sonarjs/no-dead-store": "error",
    "sonarjs/no-duplicated-branches": "error",
    "sonarjs/no-element-overwrite": "error",
    "sonarjs/no-empty-collection": "error",
    "sonarjs/no-gratuitous-expressions": "error",
    "sonarjs/no-hardcoded-passwords": "error",
    "sonarjs/no-hardcoded-secrets": "error",
    "sonarjs/no-identical-conditions": "error",
    "sonarjs/no-identical-expressions": "error",
    "sonarjs/no-identical-functions": "error",
    "sonarjs/no-ignored-exceptions": "error",
    "sonarjs/no-inverted-boolean-check": "error",
    "sonarjs/no-invariant-returns": "error",
    "sonarjs/no-nested-assignment": "error",
    "sonarjs/no-nested-switch": "error",
    "sonarjs/no-nested-template-literals": "error",
    "sonarjs/no-redundant-assignments": "error",
    "sonarjs/no-redundant-boolean": "error",
    "sonarjs/no-redundant-jump": "error",
    "sonarjs/no-same-line-conditional": "error",
    "sonarjs/no-unused-collection": "error",
    "sonarjs/no-use-of-empty-return-value": "error",
    "sonarjs/prefer-object-literal": "error",
    "sonarjs/prefer-promise-shorthand": "error",
    "sonarjs/prefer-single-boolean-return": "error",
    "sonarjs/prefer-type-guard": "error",
    "sonarjs/prefer-while": "error",
    "sonarjs/pseudo-random": "error",
    "sonarjs/redundant-type-aliases": "error",
    "sonarjs/slow-regex": "error",
    "sonarjs/updated-loop-counter": "error",
    "sonarjs/use-type-alias": "error",

    "react/react-in-jsx-scope": "off",
    "eslint/no-shadow": "off",
    "eslint/no-await-in-loop": "off",
    "eslint/no-underscore-dangle": "off",

    "typescript/no-floating-promises": "error",
    "typescript/no-misused-promises": "error",
    "typescript/no-unsafe-assignment": "error",
    "typescript/no-unsafe-argument": "error",
    "typescript/no-unsafe-call": "error",
    "typescript/no-unsafe-member-access": "error",
    "typescript/no-unsafe-return": "error",
    "typescript/only-throw-error": "error",
    "typescript/switch-exhaustiveness-check": "error",
    "typescript/await-thenable": "error",
    "typescript/no-for-in-array": "error",
    "typescript/require-await": "error",

    "import/no-cycle": "error",

    "typescript/no-explicit-any": "error",
    "typescript/no-non-null-assertion": "error",
    "typescript/consistent-type-assertions": ["error", { assertionStyle: "never" }],
    "no-restricted-globals": ["error", { name: "Date", message: "Use luxon DateTime instead of the built-in Date." }],

    "treetime/no-vague-identifiers": "error",
    "treetime/no-exec-shell-string": "error",
    "treetime/no-module-level-mutable": "error",
    "treetime/no-versioned-names": "error",
    "treetime/no-typographic-characters": "error",
    "treetime/use-class-name-helper": "error",
    "treetime/require-io-timeout": "error",
    "treetime/no-async-array-predicate": "error",
    "treetime/require-suppression-reason": "error",
    "treetime/no-unbounded-suppression": "error",
    "treetime/callers-before-callees": "error",
    "treetime/no-side-effects-in-getters": "error",

    "anti-slop/no-array-filter-map": "error",
    "anti-slop/no-conditional-empty-object-spread": "error",
    "anti-slop/no-known-value-widening": "error",
    "anti-slop/no-module-mocking": "error",
    "anti-slop/no-object-parameters": "error",
    "anti-slop/no-reduce-accumulator-copy": "error",
    "anti-slop/no-reflect-apply": "error",
    "anti-slop/no-reflect-get": "error",
    "anti-slop/no-runtime-typeof": ["error", { allowInTypeGuards: true }],
    "anti-slop/no-shape-in-symbol-names": "error",
    "anti-slop/no-unknown-parameters": "error",
    "anti-slop/no-unknown-returns": "error",
    "anti-slop/no-unknown-type-aliases": "error",
    "anti-slop/no-unsafe-dictionary-type": "error",
    "anti-slop/require-readable-spacing": "error",
    "oxc/no-accumulating-spread": "error",

    "no-restricted-imports": ["error", { patterns: IMPORT_BOUNDARY_PATTERNS }],

    "react/rules-of-hooks": "error",
    "react/checked-requires-onchange-or-readonly": "error",
    "react/display-name": "error",
    "react/jsx-no-target-blank": "error",
    "react/jsx-no-useless-fragment": "error",
    "react/no-unescaped-entities": "error",

    "typescript/ban-types": "error",
    "typescript/no-unsafe-function-type": "error",
    "typescript/no-deprecated": "error",
    "typescript/no-confusing-void-expression": ["error", { ignoreArrowShorthand: true }],
    "typescript/no-mixed-enums": "error",
    "typescript/prefer-enum-initializers": "error",
    "typescript/prefer-includes": "error",
    "typescript/prefer-nullish-coalescing": "error",
    "typescript/prefer-promise-reject-errors": "error",
    "typescript/related-getter-setter-pairs": "error",
    "typescript/restrict-plus-operands": "error",
    "typescript/return-await": "error",
    "typescript/strict-boolean-expressions": "error",
    "typescript/strict-void-return": "error",

    "eslint/accessor-pairs": "error",
    "eslint/array-callback-return": "error",
    "eslint/eqeqeq": ["error", "always", { null: "ignore" }],
    "eslint/no-array-constructor": "error",
    "eslint/no-case-declarations": "error",
    "eslint/no-constructor-return": "error",
    "eslint/no-else-return": "error",
    "eslint/no-loop-func": "error",
    "eslint/no-new-wrappers": "error",
    "eslint/no-object-constructor": "error",
    "eslint/no-promise-executor-return": "error",
    "eslint/no-prototype-builtins": "error",
    "eslint/no-useless-return": "error",
    "eslint/radix": "error",
    "eslint/require-unicode-regexp": "error",
    "eslint/symbol-description": "error",

    "unicorn/consistent-assert": "error",
    "unicorn/consistent-empty-array-spread": "error",
    "unicorn/escape-case": "error",
    "unicorn/explicit-length-check": "error",
    "unicorn/new-for-builtins": "error",
    "unicorn/no-hex-escape": "error",
    "unicorn/no-immediate-mutation": "error",
    "unicorn/no-instanceof-array": "error",
    "unicorn/no-negated-condition": "error",
    "unicorn/no-negation-in-equality-check": "error",
    "unicorn/no-new-buffer": "error",
    "unicorn/no-object-as-default-parameter": "error",
    "unicorn/no-static-only-class": "error",
    "unicorn/no-this-assignment": "error",
    "unicorn/no-typeof-undefined": "error",
    "unicorn/no-unnecessary-array-flat-depth": "error",
    "unicorn/no-unnecessary-array-splice-count": "error",
    "unicorn/no-unnecessary-slice-end": "error",
    "unicorn/no-unreadable-iife": "error",
    "unicorn/no-useless-promise-resolve-reject": "error",
    "unicorn/no-useless-switch-case": "error",
    "unicorn/prefer-array-flat": "error",
    "unicorn/prefer-array-some": "error",
    "unicorn/prefer-at": "error",
    "unicorn/prefer-blob-reading-methods": "error",
    "unicorn/prefer-code-point": "error",
    "unicorn/prefer-import-meta-properties": "error",
    "unicorn/prefer-math-min-max": "error",
    "unicorn/prefer-math-trunc": "error",
    "unicorn/prefer-native-coercion-functions": "error",
    "unicorn/prefer-number-coercion": "error",
    "unicorn/prefer-prototype-methods": "error",
    "unicorn/prefer-regexp-test": "error",
    "unicorn/prefer-single-call": "error",
    "unicorn/prefer-string-replace-all": "error",
    "unicorn/prefer-string-slice": "error",
    "unicorn/prefer-top-level-await": "error",
    "unicorn/prefer-type-error": "error",
    "unicorn/require-number-to-fixed-digits-argument": "error",

    "oxc/branches-sharing-code": "error",
  },

  overrides: [
    ...packageBoundaryOverrides(),
    {
      files: ["packages/app-web/src/**", "packages/app-ui/src/**", "packages/app-desktop/renderer/**"],
      rules: {
        "import/no-nodejs-modules": "error",
      },
    },
    {
      files: ["packages/app-ui/src/**", "packages/app-web/src/**", "packages/app-desktop/renderer/**"],
      rules: {
        "web/tailwind-classes": "error",

        "web/use-themed-cn": "error",

        "react/purity": "error",
        "react/immutability": "error",
        "react/refs": "error",
        "react/set-state-in-effect": "error",
        "react/static-components": "error",

        "react/no-deriving-state-in-effects": "error",
        "web/no-chained-state-updates": "error",
        "web/no-event-handler-effect": "error",
        "web/no-state-in-effect-initializer": "error",

        "react/no-danger": "error",
        "react/forbid-dom-props": [
          "error",
          { forbid: [{ propName: "style", message: "Use Tailwind classes, not inline styles." }] },
        ],
        "no-restricted-properties": [
          "error",
          { property: "innerHTML", message: "Assigning innerHTML injects raw HTML. Render via React." },
          { property: "outerHTML", message: "Assigning outerHTML injects raw HTML. Render via React." },
          { property: "insertAdjacentHTML", message: "insertAdjacentHTML injects raw HTML. Render via React." },
          { object: "document", property: "write", message: "document.write injects raw HTML. Render via React." },
        ],
      },
    },
    {
      files: ["packages/app-desktop/src/**"],
      rules: {
        "import/no-nodejs-modules": "off",
        "unicorn/prefer-top-level-await": "off",
      },
    },
    {
      files: ["**/vite.config.ts", "**/*.config.ts", "**/scripts/**"],
      rules: {
        "import/no-nodejs-modules": "off",
        "treetime/no-module-level-mutable": "off",
      },
    },
    {
      files: ["**/*.test.ts", "**/*.test.tsx", "**/*.spec.ts", "**/*.spec.tsx", "**/__tests__/**"],
      rules: {
        "anti-slop/no-unknown-parameters": "off",
        "anti-slop/no-unknown-returns": "off",
        "anti-slop/no-unsafe-dictionary-type": "off",
        "import/no-nodejs-modules": "error",
        "unicorn/consistent-function-scoping": "off",
        "treetime/no-fake-success": "error",
        "treetime/no-tautological-assertion": "error",
        "treetime/no-test-resource-access": "error",
        "treetime/no-assertion-in-loop": "error",
        "treetime/no-uppercase-test-title": "error",
        "treetime/prefer-strict-equal": "error",
        "treetime/no-disabled-tests": "error",
        "treetime/no-focused-tests": "error",
        "treetime/prefer-test-over-it": "error",
        "treetime/no-leaky-mocks": "error",
        "vitest/no-restricted-matchers": [
          "error",
          {
            toBeTruthy: "Assert an exact value, not truthiness.",
            toBeFalsy: "Assert an exact value, not falsiness.",
            toHaveBeenCalled: "Assert the call arguments with toHaveBeenCalledWith.",
          },
        ],
        "vitest/no-restricted-vi-methods": [
          "error",
          {
            mock: "Module mocking is banned. Inject collaborators instead.",
            doMock: "Module mocking is banned. Inject collaborators instead.",
          },
        ],
      },
    },
    {
      files: [
        "**/*.integration.test.ts",
        "**/*.integration.test.tsx",
        "**/*.process.test.ts",
        "**/*.migration.test.ts",
        "**/*.gm.test.ts",
        "**/test_gm_*.ts",
      ],
      rules: {
        "import/no-nodejs-modules": "off",
      },
    },
    {
      files: ["dev/lints/oxlint/**"],
      rules: {
        "anti-slop/no-runtime-typeof": "off",
        "anti-slop/no-unknown-parameters": "off",
        "anti-slop/no-unknown-returns": "off",
        "import/no-nodejs-modules": "off",
        "treetime/callers-before-callees": "off",
      },
    },
    {
      files: [
        "packages/app-ui/src/App.tsx",
        "packages/app-ui/src/ui/fonts.ts",
        "packages/app-web/src/main.tsx",
        "packages/app-desktop/renderer/main.tsx",
      ],
      rules: {
        "import/no-unassigned-import": "off",
      },
    },
    {
      files: ["packages/app-ui/src/**"],
      rules: {
        "react-perf/jsx-no-jsx-as-prop": "off",
        "react-perf/jsx-no-new-object-as-prop": "off",
      },
    },
    {
      files: [
        "packages/app-contracts/src/bridge.ts",
        "packages/app-desktop/src/desktop-bridge.ts",
        "packages/app-desktop/src/main.ts",
        "packages/app-web/src/bridge-web.ts",
      ],
      rules: {
        "anti-slop/no-runtime-typeof": "off",
        "anti-slop/no-unknown-parameters": "off",
        "anti-slop/no-unknown-returns": "off",
      },
    },
    {
      files: ["test/property.test.ts"],
      rules: {
        "typescript/no-unsafe-call": "off",
        "typescript/no-unsafe-member-access": "off",
        "vitest/no-standalone-expect": "off",
      },
    },
    {
      files: ["dev/lints/oxlint/rules/tailwind-classes.ts"],
      rules: {
        "unicorn/no-array-sort": "off",
      },
    },
  ],
});
