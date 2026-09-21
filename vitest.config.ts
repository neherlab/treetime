import { defineConfig } from "vitest/config";

import { DETERMINISTIC_SEED } from "./test/seed";

export default defineConfig({
  test: {
    include: ["test/**/*.test.ts", "packages/*/src/**/*.test.{ts,tsx}"],
    setupFiles: ["./vitest.setup.ts"],
    environment: "node",
    passWithNoTests: false,
    sequence: {
      shuffle: true,
      seed: DETERMINISTIC_SEED,
    },
  },
});
