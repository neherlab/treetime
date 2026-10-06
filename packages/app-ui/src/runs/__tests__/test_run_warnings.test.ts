import { describe, expect, test } from "vitest";

import { showNamesLabel } from "../RunWarnings";

describe("run warnings", () => {
  test.each([
    [1, "Show the name"],
    [2, "Show the 2 names"],
    [1200, "Show the 1200 names"],
  ])("the toggle of %i names reads %s", (count, expected) => {
    expect(showNamesLabel(count)).toBe(expected);
  });
});
