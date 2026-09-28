import { describe, expect, test } from "vitest";

import { wordMatcher } from "../text";

describe("text", () => {
  test("a query matches text holding every word, in any order and case", () => {
    const matches = wordMatcher("  Clock RATE ");

    expect([matches("rate of the clock"), matches("clock filter"), wordMatcher("")("anything")]).toStrictEqual([
      true,
      false,
      true,
    ]);
  });
});
