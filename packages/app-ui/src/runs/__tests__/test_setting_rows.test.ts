import { shownSettings, type ShownSpec } from "../settingRows";

const SPECS: ShownSpec[] = [
  { key: "tree", role: "input" },
  { key: "keep_tips_file", role: "input" },
  { key: "clock_rate", role: "setting" },
  { key: "coalescent", role: "setting" },
  { key: "outdir", role: "output" },
];

const CONFIG = { tree: "inputs/tree.nwk", clock_rate: 0.001, outdir: "out" };

const CHANGED = new Set(["clock_rate"]);

describe("settings of a run", () => {
  test("changed only lists the inputs the run read and the settings it changed", () => {
    expect(shownSettings(SPECS, CONFIG, CHANGED, true).map((spec) => spec.key)).toStrictEqual(["tree", "clock_rate"]);
  });

  test("all settings list every input and setting, never the outputs", () => {
    expect(shownSettings(SPECS, CONFIG, CHANGED, false).map((spec) => spec.key)).toStrictEqual([
      "tree",
      "keep_tips_file",
      "clock_rate",
      "coalescent",
    ]);
  });
});
