import type { RunWarning } from "@neherlab/app-contracts";

import { showNamesLabel, warningsOf, warningSummary } from "../RunWarnings";

const OUTSIDE: RunWarning = {
  kind: "dates-outside-input-dates",
  message: "3 tips dated outside their input date range.",
  names: ["A", "B", "C"],
};

const DUPLICATES: RunWarning = {
  kind: "duplicate-sequence-names",
  message: "2 sequence names occur more than once",
  names: ["D", "E"],
};

describe("run warnings", () => {
  test.each([
    [1, "Show the name"],
    [2, "Show the 2 names"],
    [1200, "Show the 1200 names"],
  ])("the toggle of %i names reads %s", (count, expected) => {
    expect(showNamesLabel(count)).toBe(expected);
  });

  test("the summary line gives the count and the first cause", () => {
    expect(warningSummary([OUTSIDE, DUPLICATES])).toBe("2 warnings: 3 tips dated outside their input date range");
  });

  test("the summary line of one warning uses the singular", () => {
    expect(warningSummary([DUPLICATES])).toBe("1 warning: 2 sequence names occur more than once");
  });

  test("a run without warnings has no summary line", () => {
    expect(warningSummary([])).toBeUndefined();
  });

  test("the samples estimate shows only the warnings about samples", () => {
    expect(warningsOf([DUPLICATES, OUTSIDE], "samples")).toStrictEqual([OUTSIDE]);
  });
});
