import { describe, expect, test } from "vitest";

import { shouldRestart } from "../backend-process";

describe("backend_process restarts", () => {
  test("the first crash restarts the back end", () => {
    expect(shouldRestart([1_000], 1_000)).toBe(true);
  });

  test("four crashes within a minute still restart the back end", () => {
    expect(shouldRestart([0, 10_000, 20_000, 30_000], 30_000)).toBe(true);
  });

  test("the fifth crash within a minute stops the restarts", () => {
    expect(shouldRestart([0, 10_000, 20_000, 30_000, 40_000], 40_000)).toBe(false);
  });

  test("crashes older than a minute do not count", () => {
    expect(shouldRestart([0, 1_000, 2_000, 3_000, 70_000], 70_000)).toBe(true);
  });
});
