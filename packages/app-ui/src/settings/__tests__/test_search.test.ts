import { describe, expect, test } from "vitest";

import { defaultConfig } from "../config";
import { setAt } from "../json";
import { commandSettings } from "../schema";
import { matchingSpecs } from "../search";

const specs = commandSettings("clock").specs;

describe("setting search", () => {
  test("matches the flag of a nested setting", () => {
    expect(
      matchingSpecs(specs, defaultConfig(specs), "--branch-split-grid", false).map((spec) => spec.key),
    ).toStrictEqual(["branch_split.n_points"]);
  });

  test("changed only keeps the changed settings", () => {
    const config = setAt(defaultConfig(specs), ["keep_root"], true);

    expect(matchingSpecs(specs, config, "", true).map((spec) => spec.key)).toStrictEqual(["keep_root"]);
  });

  test("every word must match", () => {
    expect(matchingSpecs(specs, defaultConfig(specs), "keep root nonsense", false)).toStrictEqual([]);
  });
});
