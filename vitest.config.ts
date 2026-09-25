import { defineConfig } from "vitest/config";

const DETERMINISTIC_SEED = 20260919;

export default defineConfig({
  test: {
    include: ["packages/*/src/**/*.test.{ts,tsx}", "packages/*/build/**/*.test.ts"],
    setupFiles: ["./vitest.setup.ts"],
    environment: "node",
    passWithNoTests: false,
    provide: { seed: DETERMINISTIC_SEED },
    sequence: {
      shuffle: true,
      seed: DETERMINISTIC_SEED,
    },
  },
});
