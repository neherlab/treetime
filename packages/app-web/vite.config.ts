import { readFileSync } from "node:fs";
import { homedir } from "node:os";
import { join } from "node:path";

import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { RouteStore, formatUrl, parseHostname } from "portless";
import { createLogger, defineConfig, type Plugin } from "vite";

const apiPort = process.env["TREETIME_API_PORT"] ?? "3100";

const webPort = Number(process.env["TREETIME_WEB_PORT"] ?? "5173");

const logger = createLogger("info", { allowClearScreen: false });

const originalInfo = logger.info.bind(logger);

logger.info = (msg, options) => {
  if (msg.includes("hmr") || msg.includes("page reload")) return;
  originalInfo(msg, options);
};

export default defineConfig({
  plugins: [tailwindcss(), react(), portless()],
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
