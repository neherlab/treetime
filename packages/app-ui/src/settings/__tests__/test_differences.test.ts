import { describe, expect, test } from "vitest";

import { COMMAND_SETTINGS } from "../catalog";
import { defaultConfig } from "../config";
import { settingDifferences } from "../differences";
import { setAt } from "../json";

const SPECS = COMMAND_SETTINGS.timetree.specs;

const BASE = setAt(setAt(defaultConfig(SPECS), ["tree"], "/runs/a/inputs/tree.nwk"), ["output_all"], "/runs/a/out");

const INPUT = { setting: "tree", path: "/runs/a/inputs/tree.nwk", size: 10, sha256: "abc" };

describe("settings that differ between two runs", () => {
  test("identical runs have no difference, output folders excluded", () => {
    const other = setAt(BASE, ["output_all"], "/runs/b/out");

    expect(
      settingDifferences(SPECS, { config: BASE, inputs: [INPUT] }, { config: other, inputs: [INPUT] }),
    ).toStrictEqual([]);
  });

  test("a changed setting is listed with both values, nested keys included", () => {
    const other = setAt(setAt(BASE, ["clock_filter"], 2), ["coalescent_skyline"], true);

    const differences = settingDifferences(
      SPECS,
      { config: BASE, inputs: [INPUT] },
      { config: other, inputs: [INPUT] },
    );

    expect(
      differences
        .map((difference) => [difference.spec.key, difference.first, difference.second] as const)
        .toSorted((left, right) => left[0].localeCompare(right[0])),
    ).toStrictEqual([
      ["clock_filter", 3, 2],
      ["coalescent_skyline", false, true],
    ]);
  });

  test("the same input file at another path is listed as the same content", () => {
    const moved = { ...INPUT, path: "/runs/b/inputs/tree.nwk" };

    const [difference] = settingDifferences(
      SPECS,
      { config: BASE, inputs: [INPUT] },
      { config: BASE, inputs: [moved] },
    );

    expect(difference).toMatchObject({ kind: "input", first: [INPUT.path], second: [moved.path], sameContent: true });
  });

  test("an input with other content at the same path is listed as different content", () => {
    const changed = { ...INPUT, sha256: "def" };

    const [difference] = settingDifferences(
      SPECS,
      { config: BASE, inputs: [INPUT] },
      { config: BASE, inputs: [changed] },
    );

    expect(difference).toMatchObject({ kind: "input", sameContent: false });
  });
});
