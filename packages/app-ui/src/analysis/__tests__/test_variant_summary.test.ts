import type { CheckConfigResponse, SettingDifference } from "@neherlab/app-contracts";

import { commandSettings } from "../../settings/catalog";
import { type VariantInput, variantSummary } from "../variantSummary";

const SPECS = commandSettings("timetree").settings;

const SOURCE = { title: "ebola/100, time tree", command: "timetree" } as const;

const COALESCENT = setting("coalescent", "Coalescent prior", "analysis");

const CONFIDENCE = setting("confidence", "Confidence intervals", "analysis");

const OUTPUT_PATH = setting("output_tree_nwk", "Output tree Newick", "output");

describe("variant summary", () => {
  test("an unedited variant says that no setting changed yet", () => {
    expect(variantSummary(input({ check: valid([]) }))).toBe(
      'Variant of "ebola/100, time tree". No settings changed yet.',
    );
  });

  test("the summary grows with each changed setting", () => {
    const summaries = [[COALESCENT], [COALESCENT, CONFIDENCE]].map((differences) =>
      variantSummary(input({ check: valid(differences) })),
    );

    expect(summaries).toStrictEqual([
      'Variant of "ebola/100, time tree". Changed: coalescent prior.',
      'Variant of "ebola/100, time tree". Changed: coalescent prior, confidence intervals.',
    ]);
  });

  test("settings that the app sets for every run are no change", () => {
    expect(variantSummary(input({ check: valid([OUTPUT_PATH]) }))).toBe(
      'Variant of "ebola/100, time tree". No settings changed yet.',
    );
  });

  test("another analysis is a change", () => {
    expect(variantSummary(input({ check: { ...valid([]), command: "clock" } }))).toBe(
      'Variant of "ebola/100, time tree". Changed: analysis.',
    );
  });

  test("on the page of its source the stale bar lists the changes, so the summary names only the source", () => {
    expect(
      [valid([]), valid([COALESCENT])].map((check) => variantSummary(input({ shownRunId: "source", check }))),
    ).toStrictEqual(["Variant of this run. No settings changed yet.", "Variant of this run"]);
  });

  test.each<[string, Partial<VariantInput>]>([
    ["the check is pending", { check: undefined }],
    ["the settings are invalid", { check: invalid() }],
    ["the source is unknown", { source: undefined, check: valid([COALESCENT]) }],
    ["the page shows another run", { shownRunId: "other", check: valid([COALESCENT]) }],
  ])("without a comparison with the source (%s) the summary names only the source", (_case, changes) => {
    expect(variantSummary(input(changes))).toMatch(/^Variant of "[^"]+"$/u);
  });

  test("a draft that starts from no run and the read-only options of a run are no variant", () => {
    expect([
      variantSummary(input({ fromRunId: undefined, check: valid([COALESCENT]) })),
      variantSummary(input({ viewing: true, shownRunId: "source", check: valid([COALESCENT]) })),
    ]).toStrictEqual([undefined, undefined]);
  });
});

function input(changes: Partial<VariantInput>): VariantInput {
  return {
    fromRunId: "source",
    shownRunId: undefined,
    viewing: false,
    source: SOURCE,
    check: undefined,
    specs: SPECS,
    ...changes,
  };
}

function setting(key: string, label: string, purpose: SettingDifference["purpose"]): SettingDifference {
  return { kind: "setting", key, label, group: label, purpose, second: true };
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
