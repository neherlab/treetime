import { describe, expect, test } from "vitest";

import { fromJsonFloat, nonFiniteLabel, parseNumber, parseOptionalNumber } from "../numbers";

describe("non-finite numbers", () => {
  test("the JSON strings of infinities and NaN become numbers", () => {
    expect([fromJsonFloat("inf"), fromJsonFloat("-inf"), fromJsonFloat("nan"), fromJsonFloat(1.5)]).toStrictEqual([
      Number.POSITIVE_INFINITY,
      Number.NEGATIVE_INFINITY,
      Number.NaN,
      1.5,
    ]);
  });

  test("text spellings of infinities and NaN are read in any case", () => {
    expect(["inf", "-inf", "NaN", "Infinity", "-1e-3"].map(parseNumber)).toStrictEqual([
      Number.POSITIVE_INFINITY,
      Number.NEGATIVE_INFINITY,
      Number.NaN,
      Number.POSITIVE_INFINITY,
      -0.001,
    ]);
  });

  test("an empty cell is absent and text that is not a number is an error", () => {
    expect(parseOptionalNumber("")).toBeUndefined();
    expect(() => parseNumber("abc")).toThrow('"abc" is not a number');
  });

  test("non-finite values have a readable label", () => {
    expect([Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NaN].map(nonFiniteLabel)).toStrictEqual([
      "+infinity",
      "-infinity",
      "not a number",
    ]);
  });
});
