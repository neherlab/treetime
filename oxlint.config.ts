import { defineConfig } from "oxlint";

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

  jsPlugins: ["./dev/lints/oxlint/index.ts", "./dev/lints/oxlint/web.ts"],

  categories: {
    correctness: "warn",
    suspicious: "warn",
    perf: "warn",
  },

  options: {
    typeAware: true,
  },

  rules: {
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
        "treetime/no-always-true-assertion": "error",
        "treetime/no-disabled-tests": "error",
        "treetime/no-focused-tests": "error",
        "treetime/prefer-test-over-it": "error",
        "treetime/no-uppercase-test-title": "error",
        "treetime/no-module-mocks": "error",
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
