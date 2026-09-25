import { describe, expect, test } from "vitest";

import { APP_COMMANDS } from "../commands";
import { groupedSpecs } from "../groups";
import { commandSettings, type SettingSpec } from "../schema";

function specFor(key: string, pathRole: SettingSpec["pathRole"]): SettingSpec {
  return {
    key,
    path: key.split("."),
    kind: "text",
    nullable: false,
    options: [],
    itemKind: "string",
    defaultValue: "",
    flag: `--${key}`,
    numArgs: [1, 1],
    valueDelimiter: null,
    cliValues: {},
    pathRole,
    minimum: null,
    help: "",
    more: "",
  };
}

function groupKeys(specs: readonly SettingSpec[]) {
  return groupedSpecs(specs).map(([group, members]) => [group, members.map((spec) => spec.key)]);
}

describe("setting groups", () => {
  test.each(APP_COMMANDS)("every %s setting has a group other than other", (command) => {
    const other = groupedSpecs(commandSettings(command).specs).find(([group]) => group === "Other");

    expect(other?.[1].map((spec) => spec.key) ?? []).toStrictEqual([]);
  });

  test("a setting missing from the table falls into other", () => {
    expect(groupKeys([specFor("brand_new_setting", null), specFor("clock_rate", null)])).toStrictEqual([
      ["Molecular clock", ["clock_rate"]],
      ["Other", ["brand_new_setting"]],
    ]);
  });

  test("output paths group under outputs", () => {
    expect(groupKeys([specFor("output_tree_nwk", "output")])).toStrictEqual([["Outputs", ["output_tree_nwk"]]]);
  });

  test("nested settings follow their parent key", () => {
    expect(groupKeys([specFor("branch_split.n_points", null)])).toStrictEqual([["Rooting", ["branch_split.n_points"]]]);
  });

  test("groups keep the table order and the settings the table order within a group", () => {
    expect(
      groupKeys([specFor("seed", null), specFor("clock_std_dev", null), specFor("clock_rate", null)]),
    ).toStrictEqual([
      ["Molecular clock", ["clock_rate", "clock_std_dev"]],
      ["Reproducibility", ["seed"]],
    ]);
  });
});
