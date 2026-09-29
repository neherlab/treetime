import { describe, expect, test } from "vitest";

import { formatList, parseList } from "../lists";

describe("list settings", () => {
  test("values are separated by spaces and a tab is written as an escape", () => {
    expect(formatList([",", "\t", ";"])).toStrictEqual(", \\t ;");
  });

  test("an escaped tab reads back as a tab", () => {
    expect(parseList(", \\t ;", "string")).toStrictEqual([",", "\t", ";"]);
  });

  test("number lists read numbers and extra spaces are ignored", () => {
    expect(parseList("  1.0   0.5 ", "number")).toStrictEqual([1, 0.5]);
  });

  test("an empty text is an empty list", () => {
    expect(parseList("  ", "string")).toStrictEqual([]);
  });

  test("a word in a number list is not a number, so the form reports it", () => {
    expect(parseList("1 x", "number")).toStrictEqual([1, "x"]);
  });
});
