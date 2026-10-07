import { readFileSync } from "node:fs";
import { join } from "node:path";

import { projectConfig } from "./dev/lints/oxlint/config.ts";
import { openApiEnums } from "./dev/lints/oxlint/rules/no-contract-enum-copy.ts";

const WEB_SCOPES = ["packages/app-ui/src/**", "packages/app-web/src/**", "packages/app-desktop/renderer/**"];

export default projectConfig({
  root: import.meta.dirname,

  ignorePatterns: ["dist", "dist-electron", ".turbo", ".build", "coverage", "packages/app-contracts/src/generated"],

  tailwind: { cwd: "packages/app-web", entryPoint: "src/index.css" },

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
    { files: ["packages/app-desktop/**"], dir: "app-desktop", web: false, allow: [] },
    { files: ["packages/app-desktop/renderer/**"], dir: "app-desktop", web: true, allow: [] },
    { files: ["packages/app-ui/src/**"], dir: "app-ui", web: true, allow: [] },
    { files: ["packages/app-web/src/**"], dir: "app-web", web: true, allow: [] },
  ],

  restrictionAllowances: [
    {
      files: [
        "packages/app-ui/src/host.ts",
        "packages/app-ui/src/runs/RootToTipPlot.tsx",
        "packages/app-ui/src/settings/inputs.ts",
      ],
      dir: "app-ui",
      web: true,
      allow: ["zod"],
    },
    { files: ["packages/app-ui/src/preferences/storage.ts"], dir: "app-ui", web: true, allow: ["zod", "jsonParse"] },
    { files: ["packages/app-desktop/src/napi-error.ts"], dir: "app-desktop", web: false, allow: ["jsonParse"] },
    { files: ["packages/app-desktop/src/ipc-main.ts"], dir: "app-desktop", web: false, allow: ["ipcMain"] },
    { files: ["packages/app-desktop/src/ipc-renderer.ts"], dir: "app-desktop", web: false, allow: ["ipcRenderer"] },
  ],

  contracts: {
    package: "@neherlab/app-contracts",
    enums: openApiEnums(readFileSync(join(import.meta.dirname, "packages/app-contracts/openapi.json"), "utf8")),
    files: ["packages/app-ui/src/**", "packages/app-web/src/**", "packages/app-desktop/**"],
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
      files: WEB_SCOPES,
      rules: {
        "custom/use-themed-cn": ["error", { module: "app-ui/src/ui/cn.ts" }],
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
  ],
});
