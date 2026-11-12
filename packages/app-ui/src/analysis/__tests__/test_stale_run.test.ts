import type { SettingDifference } from "@neherlab/app-contracts";

import { staleRun, staleText, type StaleRunInput } from "../staleRun";
import { invalidCheck, validCheck } from "./checkResponses";

const COALESCENT: SettingDifference = {
  kind: "setting",
  key: "coalescent",
  label: "Coalescent prior",
  group: "Coalescent prior",
  purpose: "analysis",
  second: "skyline",
};

describe("stale run", () => {
  test("a draft with the settings of the shown run is not stale", () => {
    expect(staleRun(input({ check: validCheck([]) }))).toStrictEqual({ stale: false });
  });

  test("a changed setting makes the shown run stale and names the setting", () => {
    expect(staleRun(input({ check: validCheck([COALESCENT]) }))).toStrictEqual({
      stale: true,
      changes: ["Coalescent prior"],
    });
  });

  test("another analysis makes the shown run stale", () => {
    expect(staleRun(input({ check: { ...validCheck([]), command: "clock" } }))).toStrictEqual({
      stale: true,
      changes: ["Analysis"],
    });
  });

  test.each([
    [true, { stale: true, changes: [] }],
    [false, { stale: false }],
  ])("a rejected draft is stale while edited (%p)", (edited, expected) => {
    expect(staleRun(input({ edited, check: invalidCheck() }))).toStrictEqual(expected);
  });

  test("without a shown run nothing is stale", () => {
    expect(staleRun(input({ shownRunId: undefined, check: validCheck([COALESCENT]) }))).toStrictEqual({ stale: false });
  });
});

describe("stale text", () => {
  test.each([
    [
      ["Coalescent prior", "Confidence intervals"],
      "Settings changed since this run: coalescent prior, confidence intervals. Run to update.",
    ],
    [["GTR model"], "Settings changed since this run: GTR model. Run to update."],
    [[], "Settings changed since this run. Run to update."],
  ])("names the changes %p", (changes, expected) => {
    expect(staleText(changes)).toBe(expected);
  });
});

function input(changes: Partial<StaleRunInput>): StaleRunInput {
  return { shownRunId: "r1", runCommand: "timetree", edited: true, check: undefined, ...changes };
}
