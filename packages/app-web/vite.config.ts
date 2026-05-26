import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { createLogger, defineConfig } from "vite";

const serverPort = process.env.PORT ?? "3100";

const logger = createLogger("info", { allowClearScreen: false });
const originalInfo = logger.info.bind(logger);
logger.info = (msg, options) => {
  if (msg.includes("hmr") || msg.includes("page reload")) return;
  originalInfo(msg, options);
};

export default defineConfig({
  envPrefix: ["VITE_", "TREETIME_"],
  plugins: [tailwindcss(), react()],
  customLogger: logger,
  clearScreen: false,
  server: {
    strictPort: true,
    proxy: {
      "/api": `http://localhost:${serverPort}`,
    },
  },
});
