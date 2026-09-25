import { existsSync } from "node:fs";
import { resolve } from "node:path";

import { auspice } from "@neherlab/app-ui/build/auspice-vite";
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { defineConfig, type Plugin } from "vite";
import electron from "vite-plugin-electron/simple";

import { contentSecurityPolicy } from "./build/content-security-policy";

const projectRoot = resolve(__dirname, "../..");

const napiNode = resolve(__dirname, "../app-napi/app-napi.linux-x64-gnu.node");

process.env["ELECTRON_OVERRIDE_DIST_PATH"] ??= resolve(projectRoot, "node_modules/electron/dist");

if (existsSync(napiNode)) {
  process.env["LD_PRELOAD"] = [process.env["LD_PRELOAD"], napiNode].filter(Boolean).join(":");
}

process.env["TREETIME_PROJECT_ROOT"] ??= projectRoot;

const electronArgs = process.env["ELECTRON_DISABLE_SANDBOX"] === "1" ? ["--no-sandbox"] : [];

function contentSecurityPolicyMeta(): Plugin {
  let devServer = false;

  return {
    name: "treetime-content-security-policy",
    configResolved(config) {
      devServer = config.command === "serve";
    },
    transformIndexHtml() {
      return [
        {
          tag: "meta",
          attrs: { "http-equiv": "Content-Security-Policy", content: contentSecurityPolicy(devServer) },
          injectTo: "head-prepend",
        },
      ];
    },
  };
}

export default defineConfig({
  root: "renderer",
  plugins: [
    electron({
      main: {
        entry: {
          main: resolve(__dirname, "src/main.ts"),
          backend: resolve(__dirname, "src/backend.ts"),
        },
        vite: {
          build: {
            outDir: resolve(__dirname, "dist-electron"),
            sourcemap: true,
            rolldownOptions: {
              external: ["@neherlab/app-napi"],
            },
          },
        },
        async onstart(args) {
          await args.startup([__dirname, ...electronArgs, "--enable-logging"]);
        },
      },
      preload: {
        input: resolve(__dirname, "src/preload.ts"),
        vite: {
          build: {
            outDir: resolve(__dirname, "dist-electron"),
            sourcemap: true,
          },
        },
      },
    }),
    contentSecurityPolicyMeta(),
    auspice(),
    tailwindcss(),
    react(),
  ],
  clearScreen: false,
  build: {
    outDir: resolve(__dirname, "dist"),
    emptyOutDir: true,
  },
});
