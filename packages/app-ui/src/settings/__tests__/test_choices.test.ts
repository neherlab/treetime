import type { SettingChoice } from "@neherlab/app-contracts";

import { choiceRow } from "../choices";

const CLOCK_RATE: SettingChoice = {
  choice: "clock-rate",
  keys: ["clock_rate", "clock_std_dev"],
  options: [
    { option: "estimate", patch: [{ key: "clock_rate" }, { key: "clock_std_dev" }], settings: [] },
    {
      option: "fixed",
      patch: [{ key: "clock_rate", value: 0.001 }, { key: "clock_std_dev" }],
      settings: ["clock_rate", "clock_std_dev"],
    },
  ],
};

const BEFORE = {};

const AFTER = { clock_rate: 0.001 };

describe("choice row", () => {
  test("lists the options and shows the settings of the reported option", () => {
    const row = choiceRow(CLOCK_RATE, {
      reported: [{ choice: "clock-rate", option: "fixed" }],
      picked: undefined,
      checked: AFTER,
      config: AFTER,
    });

    expect(row).toStrictEqual({
      options: ["estimate", "fixed"],
      selected: "fixed",
      shown: ["clock_rate", "clock_std_dev"],
    });
  });

  test("a picked option stays selected until the check answers for the new config", () => {
    const row = choiceRow(CLOCK_RATE, {
      reported: [{ choice: "clock-rate", option: "estimate" }],
      picked: "fixed",
      checked: BEFORE,
      config: AFTER,
    });

    expect([row.selected, row.shown]).toStrictEqual(["fixed", ["clock_rate", "clock_std_dev"]]);
  });

  test("the reported option wins once the check answers for the current config", () => {
    const row = choiceRow(CLOCK_RATE, {
      reported: [{ choice: "clock-rate", option: "estimate" }],
      picked: "fixed",
      checked: AFTER,
      config: AFTER,
    });

    expect([row.selected, row.shown]).toStrictEqual(["estimate", []]);
  });

  test("a pending check keeps the picked option", () => {
    const row = choiceRow(CLOCK_RATE, { reported: undefined, picked: "fixed", checked: undefined, config: AFTER });

    expect(row.selected).toBe("fixed");
  });

  test("without a report or a pick the first option is selected", () => {
    const row = choiceRow(CLOCK_RATE, { reported: [], picked: undefined, checked: undefined, config: BEFORE });

    expect([row.selected, row.shown]).toStrictEqual(["estimate", []]);
  });

  test("a report for another choice does not select an option of this one", () => {
    const row = choiceRow(CLOCK_RATE, {
      reported: [{ choice: "root", option: "keep" }],
      picked: undefined,
      checked: BEFORE,
      config: BEFORE,
    });

    expect(row.selected).toBe("estimate");
  });
});
