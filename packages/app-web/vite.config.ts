import { readFileSync } from "node:fs";
import { homedir } from "node:os";
import { join } from "node:path";
import zlib from "node:zlib";

import { auspice } from "@neherlab/app-ui/build/auspice-vite";
import { contentSecurityPolicyMeta } from "@neherlab/app-ui/build/content-security-policy-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { RouteStore, formatUrl, parseHostname } from "portless";
import { createLogger, defineConfig, type Plugin } from "vite";
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
    icons(),
    tailwindcss(),
    react(),
    portless(),
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

function portless(): Plugin {
  return {
    name: "treetime-portless",
    apply: "serve",
    configureServer(server) {
      const mode = process.env["TREETIME_PORTLESS"] ?? "allowed";

      if (mode === "off") return;

      const stateDir = process.env["PORTLESS_STATE_DIR"] ?? join(homedir(), ".portless");
      const store = new RouteStore(stateDir);
      const hostname = parseHostname(process.env["TREETIME_HOSTNAME"] ?? "treetime");

      server.httpServer?.once("listening", () => {
        let proxyPort: number;

        try {
          proxyPort = Number(readFileSync(store.portFilePath, "utf8").trim());
        } catch (error) {
          if (mode === "required") throw new Error(`portless proxy is not set up in ${stateDir}`, { cause: error });

          return;
        }

        store.addRoute(hostname, webPort, 0);
        server.config.logger.info(`  portless: ${formatUrl(hostname, proxyPort, true)}`);
      });
      process.once("exit", () => {
        store.removeRoute(hostname, 0);
      });
    },
  };
}
