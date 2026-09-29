import { describe, expect, test } from "vitest";

import { matchingPaletteItems, paletteGroups, paletteItem } from "../paletteItems";

const NOOP = () => undefined;

const ITEMS = [
  paletteItem("Example", "example-flu", "Influenza H3N2", "flu/h3n2/20.yaml", NOOP),
  paletteItem("Run", "run-1", "Ebola time tree", "Time tree  --clock-rate", NOOP),
  paletteItem("Action", "action-new", "New analysis", "", NOOP),
  paletteItem("Run", "run-2", "Zika clock", "Clock", NOOP),
];

describe("command palette items", () => {
  test("groups follow the order actions, runs, settings, examples, and keep the item order within a group", () => {
    expect(paletteGroups(ITEMS).map((group) => [group.kind, group.items.map((item) => item.id)])).toStrictEqual([
      ["Action", ["action-new"]],
      ["Run", ["run-1", "run-2"]],
      ["Example", ["example-flu"]],
    ]);
  });

  test("no group is listed for a kind without items", () => {
    expect(paletteGroups([]).length).toBe(0);
  });

  test("a query matches the items whose kind, title and description hold every word, in any order and case", () => {
    expect(
      ["clock-rate EBOLA", "run tree", "ebola zika", "clock"].map((query) =>
        matchingPaletteItems(ITEMS, query).map((item) => item.id),
      ),
    ).toStrictEqual([["run-1"], ["run-1"], [], ["run-1", "run-2"]]);
  });

  test("an item keeps the element to focus once the palette closes", () => {
    expect(paletteItem("Setting", "setting-a", "A", "", NOOP, "field-a").focusId).toBe("field-a");
  });
});
