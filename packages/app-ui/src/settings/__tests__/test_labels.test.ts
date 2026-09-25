import { describe, expect, test } from "vitest";

import { settingLabel } from "../labels";

describe("setting labels", () => {
  test.each([
    ["clock_rate", "Clock rate"],
    ["clock_std_dev", "Clock rate std. dev."],
    ["output_tree_nwk", "Output tree Newick"],
    ["gtr_iterations", "GTR iterations"],
    ["branch_split.n_points", "Branch split: Number of points"],
  ])("%s reads as %s", (key, label) => {
    expect(settingLabel(key)).toStrictEqual(label);
  });
});
