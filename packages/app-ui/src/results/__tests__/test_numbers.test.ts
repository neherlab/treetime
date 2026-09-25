import { describe, expect, test } from "vitest";

import { fromJsonFloat, nonFiniteLabel } from "../numbers";

describe("non-finite numbers", () => {
  test("the JSON strings of infinities and NaN become numbers", () => {
    expect([fromJsonFloat("inf"), fromJsonFloat("-inf"), fromJsonFloat("nan"), fromJsonFloat(1.5)]).toStrictEqual([
      Number.POSITIVE_INFINITY,
      Number.NEGATIVE_INFINITY,
      Number.NaN,
      1.5,
    ]);
  });

  test("non-finite values have a readable label", () => {
    expect([Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NaN].map(nonFiniteLabel)).toStrictEqual([
      "+infinity",
      "-infinity",
      "not a number",
    ]);
  });
});
