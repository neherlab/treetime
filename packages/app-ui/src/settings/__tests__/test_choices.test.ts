import type { SettingChoice } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { choiceRow } from "../choices";

const COALESCENT: SettingChoice = {
  choice: "coalescent-prior",
  keys: ["coalescent", "coalescent_opt"],
  options: [
    { option: "none", patch: [{ path: ["coalescent"] }, { path: ["coalescent_opt"] }], settings: [] },
    {
      option: "fixed",
      patch: [{ path: ["coalescent"], value: 0.5 }, { path: ["coalescent_opt"] }],
      settings: ["coalescent"],
    },
    { option: "optimized", patch: [{ path: ["coalescent"] }, { path: ["coalescent_opt"], value: true }], settings: [] },
  ],
};

const BEFORE = { coalescent_opt: false };

const AFTER = { coalescent: 0.5, coalescent_opt: false };

describe("choice row", () => {
  test("lists the options and shows the settings of the reported option", () => {
    const row = choiceRow(COALESCENT, {
      reported: [{ choice: "coalescent-prior", option: "fixed" }],
      picked: undefined,
      checked: AFTER,
      config: AFTER,
    });

    expect(row).toStrictEqual({ options: ["none", "fixed", "optimized"], selected: "fixed", shown: ["coalescent"] });
  });

  test("a picked option stays selected until the check answers for the new config", () => {
    const row = choiceRow(COALESCENT, {
      reported: [{ choice: "coalescent-prior", option: "none" }],
      picked: "fixed",
      checked: BEFORE,
      config: AFTER,
    });

    expect([row.selected, row.shown]).toStrictEqual(["fixed", ["coalescent"]]);
  });

  test("the reported option wins once the check answers for the current config", () => {
    const row = choiceRow(COALESCENT, {
      reported: [{ choice: "coalescent-prior", option: "optimized" }],
      picked: "fixed",
      checked: AFTER,
      config: AFTER,
    });

    expect([row.selected, row.shown]).toStrictEqual(["optimized", []]);
  });

  test("a pending check keeps the picked option", () => {
    const row = choiceRow(COALESCENT, { reported: undefined, picked: "optimized", checked: undefined, config: AFTER });

    expect(row.selected).toBe("optimized");
  });

  test("without a report or a pick the first option is selected", () => {
    const row = choiceRow(COALESCENT, { reported: [], picked: undefined, checked: undefined, config: BEFORE });

    expect([row.selected, row.shown]).toStrictEqual(["none", []]);
  });

  test("a report for another choice does not select an option of this one", () => {
    const row = choiceRow(COALESCENT, {
      reported: [{ choice: "root", option: "keep" }],
      picked: undefined,
      checked: BEFORE,
      config: BEFORE,
    });

    expect(row.selected).toBe("none");
  });
});
