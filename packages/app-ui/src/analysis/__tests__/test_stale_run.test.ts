import type { CheckConfigResponse, SettingDifference } from "@neherlab/app-contracts";

import { commandSettings } from "../../settings/catalog";
import { staleRun, staleText, type StaleRunInput } from "../staleRun";

const SPECS = commandSettings("timetree").settings;

const COALESCENT: SettingDifference = {
  kind: "setting",
  key: "coalescent",
  label: "Coalescent prior",
  group: "Coalescent prior",
  purpose: "analysis",
  second: "skyline",
};

const OUTPUT_PATH: SettingDifference = {
  kind: "setting",
  key: "output_tree_nwk",
  label: "Output tree Newick",
  group: "Output",
  purpose: "output",
  first: "out/tree.nwk",
  second: "out/other.nwk",
};

describe("stale run", () => {
  test("a draft with the settings of the shown run is not stale", () => {
    expect(staleRun(input({ check: valid([]) }))).toStrictEqual({ stale: false });
  });

  test("a changed setting makes the shown run stale and names the setting", () => {
    expect(staleRun(input({ check: valid([COALESCENT]) }))).toStrictEqual({
      stale: true,
      changes: ["Coalescent prior"],
    });
  });

  test("a setting that the app sets for every run never counts", () => {
    expect(staleRun(input({ check: valid([OUTPUT_PATH]) }))).toStrictEqual({ stale: false });
  });

  test("another analysis makes the shown run stale", () => {
    expect(staleRun(input({ check: { ...valid([]), command: "clock" } }))).toStrictEqual({
      stale: true,
      changes: ["Analysis"],
    });
  });

  test.each([
    [true, { stale: true, changes: [] }],
    [false, { stale: false }],
  ])("a rejected draft is stale while edited (%p)", (edited, expected) => {
    expect(staleRun(input({ edited, check: invalid() }))).toStrictEqual(expected);
  });

  test("without a shown run nothing is stale", () => {
    expect(staleRun(input({ shownRunId: undefined, check: valid([COALESCENT]) }))).toStrictEqual({ stale: false });
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
  return { shownRunId: "r1", runCommand: "timetree", edited: true, check: undefined, specs: SPECS, ...changes };
}

function valid(differences: SettingDifference[]): CheckConfigResponse {
  return {
    status: "valid",
    command: "timetree",
    config: {},
    code: {
      command_line: [],
      command_line_text: "",
      yaml: [],
      yaml_text: "",
      config_file: "",
      config_command: "",
      methods: "",
      citation: { text: "", short: "", doi: "", url: "" },
    },
    checks: [],
    choices: [],
    effective: [],
    run_differences: differences,
  };
}

function invalid(): CheckConfigResponse {
  return {
    status: "invalid",
    command: "timetree",
    message: "bad",
    causes: [],
    problems: [],
    messages: [],
    checks: [],
    choices: [],
    effective: [],
  };
}
