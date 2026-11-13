import { changedSettingCount, shownSettings, type ShownSpec } from "../settingRows";

const SPECS: ShownSpec[] = [
  { key: "tree", role: "input", set_by_app: false },
  { key: "keep_tips_file", role: "input", set_by_app: false },
  { key: "clock_rate", role: "setting", set_by_app: false },
  { key: "coalescent", role: "setting", set_by_app: false },
  { key: "outdir", role: "output", set_by_app: true },
  { key: "output_selection", role: "setting", set_by_app: true },
];

const CONFIG = { tree: "inputs/tree.nwk", clock_rate: 0.001, outdir: "out", output_selection: ["tree"] };

const CHANGED = new Set(["clock_rate", "output_selection"]);

describe("settings of a run", () => {
  test("changed only lists the inputs the run read and the settings it changed", () => {
    expect(shownSettings(SPECS, CONFIG, CHANGED, true).map((spec) => spec.key)).toStrictEqual(["tree", "clock_rate"]);
  });

  test("all settings list every input and setting, never the settings the app sets", () => {
    expect(shownSettings(SPECS, CONFIG, CHANGED, false).map((spec) => spec.key)).toStrictEqual([
      "tree",
      "keep_tips_file",
      "clock_rate",
      "coalescent",
    ]);
  });

  test("the count of changed settings leaves out the settings the app sets", () => {
    expect(changedSettingCount(SPECS, CHANGED)).toBe(1);
  });
});
