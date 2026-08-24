import { defineConfig } from "oxfmt";

export default defineConfig({
  ignorePatterns: ["dist", "dist-electron", "node_modules", "bun.lock", "*.tsbuildinfo"],
  sortPackageJson: {},
});
