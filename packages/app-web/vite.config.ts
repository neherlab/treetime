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

const serverPort = process.env.PORT ?? "3100";

export default defineConfig({
<<<<<<< HEAD
  plugins: [react()],
<<<<<<< HEAD
=======
=======
  plugins: [tailwindcss(), react()],
>>>>>>> 0d74b8c0 (feat(app-ui): add UI mockup with layout shell, input panel, and results placeholders)
  customLogger: logger,
  clearScreen: false,
>>>>>>> 3eb7c125 (chore: unify dev tooling under turbo orchestration)
  server: {
    strictPort: true,
    proxy: {
      "/api": `http://localhost:${serverPort}`,
    },
  },
});
