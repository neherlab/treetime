import type { AppCommand } from "@neherlab/app-contracts";

import { runRecord } from "../../__tests__/api_server";
import { settingText } from "../settingText";

describe("setting text of a run", () => {
  test.each<[string, AppCommand, string, Record<string, string>, string]>([
    ["the label of a coalescent option", "timetree", "coalescent", { coalescent: "opt" }, "Optimized constant"],
    ["a coalescent time scale as given", "timetree", "coalescent", { coalescent: "26w" }, "26w"],
    ["the label of the reconstruction method", "ancestral", "method_anc", { method_anc: "marginal" }, "Marginal"],
  ])("shows %s", (_, command, key, config, expected) => {
    const record = { ...runRecord("r1", "ok"), command, config };

    expect(settingText(record, key, "en")).toBe(expected);
  });
});
