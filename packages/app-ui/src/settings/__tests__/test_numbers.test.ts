import { describe, expect, test } from "vitest";

import { parseNumber } from "../numbers";

describe("number text", () => {
  test.each([
    ["3", 3],
    ["-0.5", -0.5],
    [" 1e-12 ", 1e-12],
    [".5", 0.5],
    ["2.", 2],
    ["+4E3", 4000],
  ])("%j reads as a number", (text, number) => {
    expect(parseNumber(text)).toBe(number);
  });

  test.each(["-", "1e-", "0x10", "abc", "Infinity", "1e999", "1,5"])(
    "%j stays text, so the form reports it",
    (text) => {
      expect(parseNumber(text)).toBe(text);
    },
  );
});
