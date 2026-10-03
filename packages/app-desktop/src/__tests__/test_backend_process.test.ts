import { describe, expect, test } from "vitest";

import { shouldRestart, stopReason } from "../backend-process";

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

describe("backend_process stop reasons", () => {
  test("a requested restart names the new runs folder as its cause", () => {
    expect(stopReason(0, true, true)).toBe(
      "the back end restarts to open the new runs folder; the request was not answered",
    );
  });

  test("a crash that restarts names its exit code", () => {
    expect(stopReason(3, false, true)).toBe(
      "the back end stopped with exit code 3 and restarts; the request was not answered",
    );
  });

  test("a crash that ends the restarts asks the user to restart TreeTime", () => {
    expect(stopReason(3, false, false)).toBe(
      "the back end stopped with exit code 3 too often and does not restart; restart TreeTime",
    );
  });
});
