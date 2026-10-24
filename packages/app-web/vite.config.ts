import zlib from "node:zlib";

import { auspice } from "@neherlab/app-ui/build/auspice-vite";
import { contentSecurityPolicyMeta } from "@neherlab/app-ui/build/content-security-policy-vite";
import { exampleData } from "@neherlab/app-ui/build/example-data-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { createLogger, defineConfig } from "vite";
import { compression, defineAlgorithm } from "vite-plugin-compression2";

const apiPort = process.env["TREETIME_API_PORT"] ?? "3100";

const webPort = Number(process.env["TREETIME_WEB_PORT"] ?? "5173");

const COMPRESSIBLE_FILES = /\.(html|css|js|mjs|json|svg|webmanifest)$/u;

const COMPRESSION_THRESHOLD = 1024;

const logger = createLogger("info", { allowClearScreen: false });

const originalInfo = logger.info.bind(logger);

logger.info = (msg, options) => {
  if (msg.includes("hmr") || msg.includes("page reload")) return;
  originalInfo(msg, options);
};

export default defineConfig({
  plugins: [
    contentSecurityPolicyMeta(),
    auspice(),
    exampleData(),
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
  publicDir,
  customLogger: logger,
  clearScreen: false,
  server: {
    port: webPort,
    strictPort: true,
    allowedHosts: [".localhost"],
    proxy: {
      "/api": `http://127.0.0.1:${apiPort}`,
    },
  },
});
