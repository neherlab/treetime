import { resolve } from "node:path";

import { auspice } from "@neherlab/app-ui/build/auspice-vite";
import { contentSecurityPolicyMeta } from "@neherlab/app-ui/build/content-security-policy-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";
import electron from "vite-plugin-electron/simple";

const projectRoot = resolve(__dirname, "../..");

process.env["TREETIME_PROJECT_ROOT"] ??= projectRoot;

const electronArgs = process.env["ELECTRON_DISABLE_SANDBOX"] === "1" ? ["--no-sandbox"] : [];

export default defineConfig({
  root: "renderer",
  base: "/",
  publicDir,
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
    icons(),
    tailwindcss(),
    react(),
  ],
  clearScreen: false,
  build: {
    outDir: resolve(__dirname, "dist"),
    emptyOutDir: true,
  },
});
