import { readFileSync } from "node:fs";
import { join } from "node:path";

import { projectConfig } from "./dev/lints/oxlint/config.ts";
import { openApiEnums } from "./dev/lints/oxlint/rules/no-contract-enum-copy.ts";

const WEB_SCOPES = [
  "packages/app-docs/src/**",
  "packages/app-ui/src/**",
  "packages/app-web/src/**",
  "packages/app-web2/src/**",
  "packages/app-desktop/renderer/**",
];

const APP_SCOPE = ["sharedComponents", "sharedComponentsAria", "sharedTheme"] as const;

export default projectConfig({
  root: import.meta.dirname,

  ignorePatterns: [
    "dist",
    "dist-electron",
    ".build",
    "coverage",
    "packages/app-contracts/src/generated",
    "packages/app-ui/src/i18n/generated",
    "packages/app-docs/.astro",
    "packages/app-wasm/pkg",
  ],

  tailwind: { cwd: "packages/app-web", entryPoint: "src/index.css" },

  tailwindScopes: [
    {
      files: ["packages/app-docs/src/**"],
      cwd: "packages/app-docs",
      entryPoint: "src/styles/docs.css",
      knownClasses: ["not-prose"],
    },
  ],

  webScopes: WEB_SCOPES,

  restrictions: {
    zod: {
      path: {
        name: "zod",
        allowTypeImports: true,
        message:
          "Validate back-end data with the zod schemas generated into @neherlab/app-contracts. Write a zod schema only for data that never reaches Rust, in an allowed module.",
      },
    },
    ipcMain: {
      property: {
        object: "ipcMain",
        message:
          "Register main-process channels through handle() and listen() of ipc-main.ts, which check the sender and the request.",
      },
    },
    ipcRenderer: {
      property: {
        object: "ipcRenderer",
        message:
          "Reach the main process through ipc-renderer.ts, which types every channel from the host channel table.",
      },
    },
    sharedComponents: {
      path: {
        name: "@base-ui/react",
        message: "Use the app-ui components, so the docs and the app share one implementation.",
      },
    },
    sharedComponentsAria: {
      path: {
        name: "react-aria-components",
        message: "Use the app-ui components, so the docs and the app share one implementation.",
      },
    },
    sharedTheme: {
      path: {
        name: "next-themes",
        message: "Use the app-ui components, so the docs and the app share one implementation.",
      },
    },
    jsonParse: {
      property: {
        object: "JSON",
        property: "parse",
        message:
          "Parsed JSON is untyped. Receive typed data from the generated client, or parse in an allowed boundary module.",
      },
    },
  },

  restrictedScopes: [
    { files: ["packages/app-desktop/**"], dir: "app-desktop", web: false, allow: APP_SCOPE },
    { files: ["packages/app-desktop/renderer/**"], dir: "app-desktop", web: true, allow: APP_SCOPE },
    { files: ["packages/app-docs/src/**"], dir: "app-docs", web: true, allow: [] },
    { files: ["packages/app-ui/src/**"], dir: "app-ui", web: true, allow: APP_SCOPE },
    { files: ["packages/app-web/src/**"], dir: "app-web", web: true, allow: APP_SCOPE },
    { files: ["packages/app-web2/src/**"], dir: "app-web2", web: true, allow: APP_SCOPE },
  ],

  restrictionAllowances: [
    {
      files: [
        "packages/app-docs/src/content.config.ts",
        "packages/app-docs/src/lib/pagefind.ts",
        "packages/app-docs/src/lib/versions.ts",
      ],
      dir: "app-docs",
      web: true,
      allow: ["zod"],
    },
    {
      files: ["packages/app-docs/src/translations/translation-status.ts"],
      dir: "app-docs",
      web: true,
      allow: ["zod", "jsonParse"],
    },
    {
      files: [
        "packages/app-ui/src/host.ts",
        "packages/app-ui/src/port.ts",
        "packages/app-ui/src/runs/RootToTipPlot.tsx",
        "packages/app-ui/src/settings/inputs.ts",
      ],
      dir: "app-ui",
      web: true,
      allow: [...APP_SCOPE, "zod"],
    },
    {
      files: ["packages/app-ui/src/preferences/storage.ts"],
      dir: "app-ui",
      web: true,
      allow: [...APP_SCOPE, "zod", "jsonParse"],
    },
    {
      files: ["packages/app-web2/src/backend/protocol.ts", "packages/app-web2/src/job/protocol.ts"],
      dir: "app-web2",
      web: true,
      allow: [...APP_SCOPE, "zod"],
    },
    {
      files: ["packages/app-desktop/src/napi-error.ts"],
      dir: "app-desktop",
      web: false,
      allow: [...APP_SCOPE, "jsonParse"],
    },
    {
      files: ["packages/app-desktop/src/ipc-main.ts"],
      dir: "app-desktop",
      web: false,
      allow: [...APP_SCOPE, "ipcMain"],
    },
    {
      files: ["packages/app-desktop/src/ipc-renderer.ts"],
      dir: "app-desktop",
      web: false,
      allow: [...APP_SCOPE, "ipcRenderer"],
    },
  ],

  contracts: {
    package: "@neherlab/app-contracts",
    enums: openApiEnums(readFileSync(join(import.meta.dirname, "packages/app-contracts/openapi.json"), "utf8")),
    files: ["packages/app-ui/src/**", "packages/app-web/src/**", "packages/app-web2/src/**", "packages/app-desktop/**"],
  },

  overrides: [
    {
      files: [
        "packages/app-desktop/src/**",
        "packages/app-ui/src/runs/ShiftPlot.tsx",
        "packages/app-ui/src/preferences/storage.ts",
      ],
      rules: {
        "custom/no-contract-parse": "off",
      },
    },
    {
      files: ["packages/app-desktop/renderer/main.tsx"],
      rules: {
        "custom/no-declare-global": "off",
      },
    },
    {
      files: [
        "packages/app-docs/src/integrations/**",
        "packages/app-docs/src/markdown/**",
        "packages/app-docs/src/preview/**",
        "packages/app-docs/src/translations/**",
      ],
      rules: {
        "import/no-nodejs-modules": "off",
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
      files: [
        "packages/app-docs/src/**/__tests__/*.test.tsx",
        "packages/app-docs/src/layouts/DocsLayout.astro",
        "packages/app-ui/src/App.tsx",
        "packages/app-ui/src/ui/fonts.ts",
        "packages/app-web/src/main.tsx",
        "packages/app-web2/src/main.tsx",
        "packages/app-desktop/renderer/main.tsx",
      ],
      rules: {
        "import/no-unassigned-import": "off",
      },
    },
    {
      files: ["packages/app-docs/src/**", "packages/app-ui/src/**"],
      rules: {
        "react-perf/jsx-no-jsx-as-prop": "off",
        "react-perf/jsx-no-new-object-as-prop": "off",
      },
    },
  ],
});
