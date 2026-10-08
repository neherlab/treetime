import zlib from "node:zlib";

import { auspice } from "@neherlab/app-ui/build/auspice-vite";
import { contentSecurityPolicyMeta } from "@neherlab/app-ui/build/content-security-policy-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";
import { compression, defineAlgorithm } from "vite-plugin-compression2";

const devPort = Number(process.env["TREETIME_WEB2_PORT"] ?? "4620");

const previewPort = Number(process.env["TREETIME_WEB2_SERVE_PORT"] ?? "2580");

const COMPRESSIBLE_FILES = /\.(html|css|js|mjs|json|svg|webmanifest|wasm)$/u;

const COMPRESSION_THRESHOLD = 1024;

export default defineConfig({
  plugins: [
    contentSecurityPolicyMeta({ webAssembly: true }),
    auspice(),
    icons(),
    tailwindcss(),
    react(),
    compression({
      include: COMPRESSIBLE_FILES,
      threshold: COMPRESSION_THRESHOLD,
      algorithms: [
        defineAlgorithm("gzip", { level: 9 }),
        defineAlgorithm("brotliCompress", { params: { [zlib.constants.BROTLI_PARAM_QUALITY]: 9 } }),
      ],
    }),
  ],
  worker: {
    format: "es",
  },
  publicDir,
  clearScreen: false,
  server: {
    port: devPort,
    strictPort: true,
    allowedHosts: [".localhost"],
  },
  preview: {
    port: previewPort,
    strictPort: true,
  },
});
