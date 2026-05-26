<<<<<<< HEAD
=======
import { resolve } from "node:path";
import { existsSync } from "node:fs";
import tailwindcss from "@tailwindcss/vite";
>>>>>>> cfc09d41 (feat(desktop): replace manual build with vite-plugin-electron)
import react from "@vitejs/plugin-react";
import electron from "vite-plugin-electron/simple";
import { defineConfig } from "vite";

const projectRoot = resolve(__dirname, "../..");
const napiNode = resolve(__dirname, "../app-napi/app-napi.linux-x64-gnu.node");

process.env.ELECTRON_OVERRIDE_DIST_PATH ??= resolve(projectRoot, "node_modules/electron/dist");

if (existsSync(napiNode)) {
  process.env.LD_PRELOAD = [process.env.LD_PRELOAD, napiNode].filter(Boolean).join(":");
}
process.env.TREETIME_PROJECT_ROOT ??= projectRoot;

export default defineConfig({
  root: "renderer",
<<<<<<< HEAD
  plugins: [react()],
=======
  plugins: [
    electron({
      main: {
        entry: resolve(__dirname, "src/main.ts"),
        vite: {
          build: {
            outDir: resolve(__dirname, "dist-electron"),
            rollupOptions: {
              external: ["@neherlab/app-napi"],
            },
          },
        },
        onstart(args) {
          args.startup([".", "--no-sandbox", "--enable-logging", "--remote-debugging-port=9229"]);
        },
      },
      preload: {
        input: resolve(__dirname, "src/preload.ts"),
        vite: {
          build: {
            outDir: resolve(__dirname, "dist-electron"),
          },
        },
      },
    }),
    tailwindcss(),
    react(),
  ],
>>>>>>> cfc09d41 (feat(desktop): replace manual build with vite-plugin-electron)
  clearScreen: false,
  server: {
    port: 5174,
    strictPort: true,
  },
  build: {
    outDir: resolve(__dirname, "dist"),
    emptyOutDir: true,
  },
});
