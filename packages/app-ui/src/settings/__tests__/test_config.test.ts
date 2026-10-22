import type { SettingSpec } from "@neherlab/app-contracts";

import { commandSettings } from "../catalog";
import {
  carryOverConfig,
  changedSpecs,
  defaultConfig,
  effectiveValue,
  isChanged,
  normalizeConfig,
  resetValue,
} from "../config";
import { withSetting } from "../json";

const clock = commandSettings("clock").settings;

const timetree = commandSettings("timetree").settings;

function spec(specs: readonly SettingSpec[], key: string): SettingSpec {
  const found = specs.find((candidate) => candidate.key === key);

  if (found === undefined) {
    throw new Error(`no setting ${key}`);
  }

  return found;
}

describe("form config", () => {
  test("defaults fill settings from the catalog", () => {
    const config = defaultConfig(clock);

    expect({ clockFilter: config["clock_filter"], keepRoot: config["keep_root"] }).toStrictEqual({
      clockFilter: 4,
      keepRoot: false,
    });
  });

  test("a default config has no changed setting", () => {
    expect(changedSpecs(timetree, defaultConfig(timetree))).toStrictEqual([]);
  });

  test("a fresh config leaves the settings with a documented default unset", () => {
    const config = defaultConfig(timetree);
    const keys = ["model", "representation", "partition_rates", "n_skyline"];

    expect([
      keys.filter((key) => Object.hasOwn(config, key)),
      keys.map((key) => effectiveValue(config, spec(timetree, key))),
    ]).toStrictEqual([[], ["infer", "auto", "estimate", 20]]);
  });

  test("a setting set to its documented default is a change", () => {
    const config = withSetting(defaultConfig(timetree), "model", "infer");

    expect(changedSpecs(timetree, config).map((found) => found.key)).toStrictEqual(["model"]);
  });

  test("a changed value is detected by its key", () => {
    const config = withSetting(defaultConfig(clock), "clock_filter", 2);

    expect(changedSpecs(clock, config).map((found) => found.key)).toStrictEqual(["clock_filter"]);
  });

  test("a list in another order is a change", () => {
    const config = withSetting(defaultConfig(timetree), "metadata_id_columns", ["name", "strain", "accession"]);

    expect(isChanged(config, spec(timetree, "metadata_id_columns"))).toStrictEqual(true);
  });

  test("input paths are never reported as changed", () => {
    const config = withSetting(defaultConfig(timetree), "tree", "data/zika/86/tree.nwk");

    expect(changedSpecs(timetree, config)).toStrictEqual([]);
  });

  test("reset restores the default of one setting only", () => {
    const changed = withSetting(withSetting(defaultConfig(timetree), "clock_rate", 0.001), "max_iter", 5);
    const clockRate = spec(timetree, "clock_rate");
    const reset = withSetting(changed, clockRate.key, resetValue(clockRate));

    expect(changedSpecs(timetree, reset).map((found) => found.key)).toStrictEqual(["max_iter"]);
  });

  test("a fresh config, a reset, and an unset setting hold no null", () => {
    const fresh = defaultConfig(timetree);
    const clockRate = spec(timetree, "clock_rate");
    const unset = withSetting(withSetting(fresh, "clock_rate", 0.003), "clock_rate", resetValue(clockRate));

    expect({
      fresh: JSON.stringify(fresh).includes("null"),
      unset: JSON.stringify(unset).includes("null"),
      hasClockRate: Object.hasOwn(unset, "clock_rate"),
      clockRateDefault: resetValue(clockRate),
      nullable: clockRate.nullable,
    }).toStrictEqual({ fresh: false, unset: false, hasClockRate: false, clockRateDefault: undefined, nullable: true });
  });

  test("normalizing drops output paths and unknown keys", () => {
    const config = normalizeConfig(clock, {
      tree: "t.nwk",
      output_all: "/elsewhere",
      output_tree_nwk: "x.nwk",
      not_a_setting: 1,
      clock_filter: 2,
    });

    expect({
      tree: config["tree"],
      outputAll: config["output_all"],
      outputTree: config["output_tree_nwk"],
      unknown: config["not_a_setting"],
      clockFilter: config["clock_filter"],
    }).toStrictEqual({
      tree: "t.nwk",
      outputAll: undefined,
      outputTree: undefined,
      unknown: undefined,
      clockFilter: 2,
    });
  });

  test("switching the command keeps inputs and changed shared settings", () => {
    const previous = withSetting(
      withSetting(withSetting(defaultConfig(timetree), "tree", "t.nwk"), "keep_root", true),
      "confidence",
      true,
    );

    const config = carryOverConfig(clock, timetree, previous);

    expect({
      tree: config["tree"],
      keepRoot: config["keep_root"],
      confidence: config["confidence"],
      changed: changedSpecs(clock, config).map((found) => found.key),
    }).toStrictEqual({ tree: "t.nwk", keepRoot: true, confidence: undefined, changed: ["keep_root"] });
  });

  test("switching the command keeps a shared choice only when the target command offers its value", () => {
    const ancestral = commandSettings("ancestral").settings;
    const homoplasy = commandSettings("homoplasy").settings;

    const previous = withSetting(
      withSetting(defaultConfig(ancestral), "output_selection", ["nwk", "gtr"]),
      "output_nwk_style",
      ["nhx"],
    );

    const config = carryOverConfig(homoplasy, ancestral, previous);

    expect({ selection: config["output_selection"], style: config["output_nwk_style"] }).toStrictEqual({
      selection: [],
      style: ["nhx"],
    });
  });
});
