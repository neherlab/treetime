import { describe, expect, test } from "vitest";

import { sortDescriptor } from "../sortDescriptor";

describe("sort descriptor", () => {
  test("an ascending sort names its column and direction", () => {
    expect(sortDescriptor([{ id: "branches", desc: false }])).toStrictEqual({
      column: "branches",
      direction: "ascending",
    });
  });

  test("a descending sort names its column and direction", () => {
    expect(sortDescriptor([{ id: "count", desc: true }])).toStrictEqual({ column: "count", direction: "descending" });
  });

  test("an unsorted table has no descriptor", () => {
    expect(sortDescriptor([])).toBeUndefined();
  });
});
