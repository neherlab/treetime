import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

const serverPort = process.env.PORT ?? "3100";

export default defineConfig({
  plugins: [react()],
  server: {
    strictPort: true,
    proxy: {
      "/api": `http://localhost:${serverPort}`,
    },
  },
});
