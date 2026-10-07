import { readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";

import { defineConfig } from "oxlint";
import type { OxlintConfig, OxlintOverride } from "oxlint";
import * as z from "zod";

const WORKSPACE_SCOPE = "@neherlab/";

const IMPORT_BOUNDARY_PATTERNS = [
  {
    group: ["@neherlab/*/src/*", "@neherlab/*/src/**", "@neherlab/*/dist/*", "@neherlab/*/dist/**"],
    message: "Import a workspace package through its entry point, not a deep internal path.",
  },
];

const SIBLING_PATH_MESSAGE =
  "A relative path must not reach into a sibling package. Import it by its @neherlab/* name.";

const REJECTED_LIBRARY_PATTERNS = [
  {
    group: ["lucide-react", "lucide-react/*"],
    message: 'Import icons from Iconify sets through unplugin-icons, e.g. `import XIcon from "~icons/lucide/x"`.',
  },
];

const PACKAGE_GRAPH_MESSAGE =
  "This package may import only the workspace packages it declares. Add the dependency to its package.json, or route through an allowed package.";

const WEB_PROPERTIES = [
  { property: "innerHTML", message: "Assigning innerHTML injects raw HTML. Render via React." },
  { property: "outerHTML", message: "Assigning outerHTML injects raw HTML. Render via React." },
  { property: "insertAdjacentHTML", message: "insertAdjacentHTML injects raw HTML. Render via React." },
  { object: "document", property: "write", message: "document.write injects raw HTML. Render via React." },
];

const TEST_FILES = ["**/*.test.ts", "**/*.test.tsx", "**/__tests__/**"];

const packageManifestSchema = z.object({
  name: z.string().optional(),
  dependencies: z.record(z.string(), z.string()).optional(),
  devDependencies: z.record(z.string(), z.string()).optional(),
  peerDependencies: z.record(z.string(), z.string()).optional(),
});

export interface ProjectSettings<R extends string> {
  root: string;
  ignorePatterns: string[];
  tailwind: { cwd: string; entryPoint: string };
  webScopes: string[];
  restrictions: Record<R, Restriction>;
  restrictedScopes: ReadonlyArray<RestrictedFiles<R>>;
  restrictionAllowances: ReadonlyArray<RestrictedFiles<R>>;
  contracts: ContractSettings | undefined;
  overrides: OxlintOverride[];
}

export type Restriction =
  | { path: { name: string; allowTypeImports?: boolean; message: string } }
  | { property: { object: string; property?: string; message: string } };

export interface RestrictedFiles<R extends string> {
  files: string[];
  dir: string;
  web: boolean;
  allow: readonly R[];
}

export interface ContractSettings {
  package: string;
  enums: ReadonlyArray<readonly string[]>;
  files: string[];
}

export function projectConfig<R extends string>(settings: ProjectSettings<R>): OxlintConfig {
  return defineConfig({
    ignorePatterns: [
      ...settings.ignorePatterns,
      "node_modules",
      "bun.lock",
      "*.tsbuildinfo",
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
      "eslint-plugin-better-tailwindcss",
      "eslint-plugin-react-you-might-not-need-an-effect",
      "./dev/lints/oxlint/index.ts",
      "./dev/lints/oxlint-anti-slop/index.ts",
    ],

    settings: {
      "better-tailwindcss": settings.tailwind,
    },

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
        {
          allow: [
            "oxlint-",
            "@ts-expect-error",
            "\\* @internal",
            "/ <reference",
            "/usr/bin/env",
            "TODO",
            "FIXME",
            "HACK",
          ],
        },
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
      "typescript/unbound-method": "error",

      "import/no-cycle": "error",

      "typescript/no-explicit-any": "error",
      "typescript/no-non-null-assertion": "error",
      "typescript/consistent-type-assertions": ["error", { assertionStyle: "never" }],
      "no-restricted-globals": ["error", { name: "Date", message: "Use luxon DateTime instead of the built-in Date." }],

      "custom/no-vague-identifiers": "error",
      "custom/no-exec-shell-string": "error",
      "custom/no-module-level-mutable": "error",
      "custom/no-versioned-names": "error",
      "custom/no-typographic-characters": "error",
      "custom/use-class-name-helper": "error",
      "custom/require-io-timeout": "error",
      "custom/no-async-array-predicate": "error",
      "custom/require-suppression-reason": "error",
      "custom/no-unbounded-suppression": "error",
      "custom/callers-before-callees": "error",
      "custom/no-side-effects-in-getters": "error",
      "custom/no-declare-global": "error",

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

      "no-restricted-imports": ["error", { patterns: [...IMPORT_BOUNDARY_PATTERNS, ...REJECTED_LIBRARY_PATTERNS] }],

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
      ...packageBoundaryOverrides(settings.root),
      ...restrictionOverrides(settings),
      ...contractOverrides(settings.contracts),
      {
        files: settings.webScopes,
        rules: {
          "import/no-nodejs-modules": "error",

          "better-tailwindcss/no-unknown-classes": "error",
          "better-tailwindcss/no-conflicting-classes": "error",
          "better-tailwindcss/no-duplicate-classes": "error",

          "react/purity": "error",
          "react/immutability": "error",
          "react/refs": "error",
          "react/set-state-in-effect": "error",
          "react/static-components": "error",

          "react/no-deriving-state-in-effects": "error",
          "react-you-might-not-need-an-effect/no-adjust-state-on-prop-change": "error",
          "react-you-might-not-need-an-effect/no-chain-state-updates": "error",
          "react-you-might-not-need-an-effect/no-derived-state": "error",
          "react-you-might-not-need-an-effect/no-event-handler": "error",
          "react-you-might-not-need-an-effect/no-external-store-subscription": "error",
          "react-you-might-not-need-an-effect/no-initialize-state": "error",
          "react-you-might-not-need-an-effect/no-pass-data-to-parent": "error",
          "react-you-might-not-need-an-effect/no-pass-live-state-to-parent": "error",
          "react-you-might-not-need-an-effect/no-reset-all-state-on-prop-change": "error",

          "react/no-danger": "error",
          "react/forbid-dom-props": [
            "error",
            { forbid: [{ propName: "style", message: "Use Tailwind classes, not inline styles." }] },
          ],
        },
      },
      ...settings.overrides,
      {
        files: ["**/vite.config.ts", "**/*.config.ts", "**/scripts/**"],
        rules: {
          "import/no-nodejs-modules": "off",
          "custom/no-module-level-mutable": "off",
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
          "custom/no-fake-success": "error",
          "custom/no-tautological-assertion": "error",
          "custom/no-test-resource-access": "error",
          "custom/no-assertion-in-loop": "error",
          "custom/no-uppercase-test-title": "error",
          "custom/prefer-strict-equal": "error",
          "custom/no-disabled-tests": "error",
          "custom/no-focused-tests": "error",
          "custom/prefer-test-over-it": "error",
          "custom/no-leaky-mocks": "error",
          "custom/no-contract-alias": "off",
          "custom/no-contract-enum-copy": "off",
          "custom/no-contract-parse": "off",
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
          "custom/callers-before-callees": "off",
        },
      },
    ],
  });
}

export function packageGraphPatterns(root: string, dir: string): Array<{ group: string[]; message: string }> {
  const packages = workspacePackages(root);
  const entry = packages.find((candidate) => candidate.dir === dir);

  if (entry === undefined) {
    return [];
  }

  const allowed = new Set([entry.name, ...entry.dependencies]);
  const banned = packages.map((other) => other.name).filter((name) => !allowed.has(name));

  return banned.length === 0 ? [] : [{ group: banned, message: PACKAGE_GRAPH_MESSAGE }];
}

function packageBoundaryOverrides(root: string): OxlintOverride[] {
  return workspacePackages(root).map((entry) => ({
    files: [`packages/${entry.dir}/**`],
    rules: { "no-restricted-imports": importRule(root, entry.dir, []) },
  }));
}

function restrictionOverrides<R extends string>(settings: ProjectSettings<R>): OxlintOverride[] {
  const names = Object.keys(settings.restrictions).filter((key): key is R => key in settings.restrictions);

  const override = (files: string[], { dir, web, allow }: RestrictedFiles<R>): OxlintOverride => {
    const applied = names.filter((name) => !allow.includes(name)).map((name) => settings.restrictions[name]);

    return {
      files,
      rules: {
        "no-restricted-imports": importRule(settings.root, dir, applied),
        "no-restricted-properties": propertyRule(web, applied),
      },
    };
  };

  return [
    ...settings.restrictedScopes.map((scope) => override(scope.files, scope)),
    ...settings.restrictionAllowances.map((allowance) => override(allowance.files, allowance)),
    ...settings.restrictedScopes.map((scope) =>
      override(
        scope.files.flatMap((files) => TEST_FILES.map((tests) => `${files.replace(/\*\*$/u, "")}${tests}`)),
        { ...scope, allow: names },
      ),
    ),
  ];
}

function contractOverrides(contracts: ContractSettings | undefined): OxlintOverride[] {
  if (contracts === undefined) {
    return [];
  }

  return [
    {
      files: contracts.files,
      rules: {
        "custom/no-contract-alias": ["error", { package: contracts.package }],
        "custom/no-contract-enum-copy": ["error", { enums: contracts.enums }],
        "custom/no-contract-parse": ["error", { package: contracts.package }],
      },
    },
  ];
}

function importRule(root: string, dir: string, applied: readonly Restriction[]): ImportRule {
  const paths = applied.flatMap((restriction) => ("path" in restriction ? [restriction.path] : []));

  const siblings = readdirSync(join(root, "packages"))
    .filter((other) => other !== dir)
    .map((other) => other.replaceAll(/[.*+?^${}()|[\]\\]/gu, String.raw`\$&`));

  const patterns = [
    ...IMPORT_BOUNDARY_PATTERNS,
    ...REJECTED_LIBRARY_PATTERNS,
    { regex: `^(\\.\\./)+(${siblings.join("|")})(/|$)`, message: SIBLING_PATH_MESSAGE },
    ...packageGraphPatterns(root, dir),
  ];

  return ["error", { paths, patterns }];
}

function propertyRule(web: boolean, applied: readonly Restriction[]): PropertyRule {
  const properties = [
    ...(web ? WEB_PROPERTIES : []),
    ...applied.flatMap((restriction) => ("property" in restriction ? [restriction.property] : [])),
  ];

  return properties.length === 0 ? "off" : ["error", ...properties];
}

type ImportRule = NonNullable<NonNullable<OxlintOverride["rules"]>["no-restricted-imports"]>;

type PropertyRule = NonNullable<NonNullable<OxlintOverride["rules"]>["no-restricted-properties"]>;

interface WorkspacePackage {
  dir: string;
  name: string;
  dependencies: string[];
}

function workspacePackages(root: string): WorkspacePackage[] {
  const packagesDir = join(root, "packages");
  const packages: WorkspacePackage[] = [];

  for (const dir of readdirSync(packagesDir)) {
    let raw: string;

    try {
      raw = readFileSync(join(packagesDir, dir, "package.json"), "utf8");
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
