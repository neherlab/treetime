import { defineConfig } from "oxlint";

const PACKAGE_GRAPH: Record<string, readonly string[]> = {
  "app-contracts": [],
  "app-napi": [],
  "app-ui": ["app-contracts"],
  "app-web": ["app-contracts", "app-ui"],
  "app-desktop": ["app-contracts", "app-napi", "app-ui"],
};

const PACKAGE_GRAPH_MESSAGE =
  "This package may import only the workspace packages it declares. Add the dependency to package.json and PACKAGE_GRAPH, or route through an allowed package.";

function packageBoundaryOverrides() {
  const names = Object.keys(PACKAGE_GRAPH);
  return names.flatMap((name) => {
    const allowed = new Set([name, ...(PACKAGE_GRAPH[name] ?? [])]);
    const banned = names.filter((other) => !allowed.has(other)).map((other) => `@neherlab/${other}`);
    if (banned.length === 0) {
      return [];
    }
    return [
      {
        files: [`packages/${name}/**`],
        rules: {
          "no-restricted-imports": ["error", { patterns: [{ group: banned, message: PACKAGE_GRAPH_MESSAGE }] }],
        },
      },
    ];
  });
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
    perf: "warn",
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
    "typescript/require-await": "warn",

    "import/no-cycle": "error",

    "typescript/no-explicit-any": "error",
    "typescript/no-non-null-assertion": "error",
    "typescript/consistent-type-assertions": ["error", { assertionStyle: "never" }],
    "no-restricted-globals": [
      "error",
      { name: "Date", message: "Use luxon DateTime instead of the built-in Date." },
    ],

    "treetime/no-vague-identifiers": "error",
    "treetime/no-exec-shell-string": "error",
    "treetime/no-module-level-mutable": "error",
    "treetime/no-versioned-names": "error",
    "treetime/no-typographic-characters": "error",
    "treetime/use-class-name-helper": "error",
    "treetime/require-io-timeout": "error",
    "treetime/no-async-array-predicate": "error",

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

    "no-restricted-imports": [
      "error",
      {
        patterns: [
          {
            group: [
              "@neherlab/*/src/*",
              "@neherlab/*/src/**",
              "@neherlab/*/dist/*",
              "@neherlab/*/dist/**",
            ],
            message: "Import a workspace package through its entry point, not a deep internal path.",
          },
          {
            group: ["../../*", "../../**"],
            message:
              "A relative path must not reach into a sibling package. Import it by its @neherlab/* name.",
          },
        ],
      },
    ],
  },

  overrides: [
    ...packageBoundaryOverrides(),
    {
      files: [
        "packages/app-web/src/**",
        "packages/app-ui/src/**",
        "packages/app-desktop/renderer/**",
      ],
      rules: {
        "import/no-nodejs-modules": "error",
      },
    },
    {
      files: ["packages/app-ui/src/**", "packages/app-web/src/**"],
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
      },
    },
    {
      files: [
        "**/vite.config.ts",
        "**/*.config.ts",
        "**/scripts/**",
        "packages/app-napi/**",
      ],
      rules: {
        "import/no-nodejs-modules": "off",
        "treetime/no-module-level-mutable": "off",
      },
    },
    {
      files: [
        "**/*.test.ts",
        "**/*.test.tsx",
        "**/*.spec.ts",
        "**/*.spec.tsx",
        "**/__tests__/**",
      ],
      rules: {
        "import/no-nodejs-modules": "error",
        "treetime/no-fake-success": "error",
        "treetime/no-tautological-assertion": "error",
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
  ],
});
