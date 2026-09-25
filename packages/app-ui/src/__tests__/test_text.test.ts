import { describe, expect, test } from "vitest";

import { truncate, wordMatcher } from "../text";

describe("text", () => {
  test("a query matches text holding every word, in any order and case", () => {
    const matches = wordMatcher("  Clock RATE ");

    expect([matches("rate of the clock"), matches("clock filter"), wordMatcher("")("anything")]).toStrictEqual([
      true,
      false,
      true,
    ]);
  });

  test("text longer than the limit ends in an ellipsis within the limit", () => {
    expect([truncate("abcdefghij", 10), truncate("abcdefghijk", 10)]).toStrictEqual(["abcdefghij", "abcdefg..."]);
  });
});
