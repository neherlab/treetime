import { defineConfig } from "vitest/config";

import { auspice } from "./packages/app-ui/build/auspice-vite.ts";
import { icons } from "./packages/app-ui/build/icons-vite.ts";

const DETERMINISTIC_SEED = 20260919;

export default defineConfig({
  plugins: [auspice(), icons()],
  test: {
    server: { deps: { inline: ["auspice"] } },
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
