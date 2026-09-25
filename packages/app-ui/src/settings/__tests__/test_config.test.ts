import { describe, expect, test } from "vitest";

import { carryOverConfig, changedSpecs, defaultConfig, isChanged, normalizeConfig, resetValue } from "../config";
import { setAt } from "../json";
import { commandSettings, type SettingSpec } from "../schema";

const clock = commandSettings("clock").specs;

const timetree = commandSettings("timetree").specs;

function spec(specs: readonly SettingSpec[], key: string): SettingSpec {
  const found = specs.find((candidate) => candidate.key === key);

  if (found === undefined) {
    throw new Error(`no setting ${key}`);
  }

  return found;
}

describe("form config", () => {
  test("defaults fill nested settings from the schema", () => {
    const config = defaultConfig(clock);

    expect({ branchSplit: config["branch_split"], clockFilter: config["clock_filter"] }).toStrictEqual({
      branchSplit: {
        method: "grid",
        n_points: 11,
        brent_max_iters: 50,
        brent_tolerance: 1e-12,
        golden_max_iters: 50,
        golden_tolerance: 1e-12,
      },
      clockFilter: 3,
    });
  });

  test("a default config has no changed setting", () => {
    expect(changedSpecs(timetree, defaultConfig(timetree))).toStrictEqual([]);
  });

  test("a changed nested value is detected by its full key", () => {
    const config = setAt(defaultConfig(clock), ["branch_split", "n_points"], 21);

    expect(changedSpecs(clock, config).map((found) => found.key)).toStrictEqual(["branch_split.n_points"]);
  });

  test("a list in another order is a change", () => {
    const config = setAt(defaultConfig(timetree), ["metadata_id_columns"], ["name", "strain", "accession"]);

    expect(isChanged(config, spec(timetree, "metadata_id_columns"))).toStrictEqual(true);
  });

  test("input paths are never reported as changed", () => {
    const config = setAt(defaultConfig(timetree), ["tree"], "data/zika/86/tree.nwk");

    expect(changedSpecs(timetree, config)).toStrictEqual([]);
  });

  test("reset restores the default of one setting only", () => {
    const changed = setAt(setAt(defaultConfig(timetree), ["clock_rate"], 0.001), ["max_iter"], 5);
    const clockRate = spec(timetree, "clock_rate");
    const reset = setAt(changed, clockRate.path, resetValue(clockRate));

    expect(changedSpecs(timetree, reset).map((found) => found.key)).toStrictEqual(["max_iter"]);
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
    }).toStrictEqual({ tree: "t.nwk", outputAll: null, outputTree: null, unknown: undefined, clockFilter: 2 });
  });

  test("switching the command keeps inputs and changed shared settings", () => {
    const previous = setAt(
      setAt(setAt(defaultConfig(timetree), ["tree"], "t.nwk"), ["keep_root"], true),
      ["confidence"],
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
});
