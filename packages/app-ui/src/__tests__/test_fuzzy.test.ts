import { describe, expect, test } from "vitest";

import { fuzzyFilter } from "../fuzzy";

const TEXTS = ["rate of the clock", "clock filter", "ancestral reconstruction", "time tree"];

function matching(query: string): string[] {
  return fuzzyFilter(TEXTS, query, (text) => text);
}

describe("fuzzy filter", () => {
  test("a query matches text holding every word, in any order and case", () => {
    expect(matching("  Clock RATE ")).toStrictEqual(["rate of the clock"]);
  });

  test("matches keep the order of the input", () => {
    expect(matching("clock")).toStrictEqual(["rate of the clock", "clock filter"]);
  });

  test("a query without letters or digits keeps every item", () => {
    expect([matching(""), matching("  "), matching("--")]).toStrictEqual([TEXTS, TEXTS, TEXTS]);
  });

  test("a word of five or more letters matches with one deleted, substituted, or swapped letter", () => {
    expect([matching("ancestrl"), matching("ancxstral"), matching("ancsetral")]).toStrictEqual([
      ["ancestral reconstruction"],
      ["ancestral reconstruction"],
      ["ancestral reconstruction"],
    ]);
  });

  test("a word with two errors does not match", () => {
    expect(matching("ancxstrl")).toStrictEqual([]);
  });

  test("a word of one or two letters matches only exactly", () => {
    expect([matching("ee"), matching("tm")]).toStrictEqual([["time tree"], []]);
  });

  test("a word prefixed with a minus excludes the items holding it", () => {
    expect(matching("clock -rate")).toStrictEqual(["clock filter"]);
  });

  test("an empty list stays empty", () => {
    expect(fuzzyFilter([], "clock", (text: string) => text)).toStrictEqual([]);
  });
});
