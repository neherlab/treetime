import type { SettingDifference } from "@neherlab/app-contracts";

import { type VariantInput, variantSummary } from "../variantSummary";
import { invalidCheck, validCheck } from "./checkResponses";

const SOURCE = { title: "ebola/100, time tree", command: "timetree" } as const;

const COALESCENT = setting("coalescent", "Coalescent prior", "analysis");

const CONFIDENCE = setting("confidence", "Confidence intervals", "analysis");

describe("variant summary", () => {
  test("an unedited variant says that no setting changed yet", () => {
    expect(variantSummary(input({ check: validCheck([]) }))).toBe(
      'Variant of "ebola/100, time tree". No settings changed yet.',
    );
  });

  test("the summary grows with each changed setting", () => {
    const summaries = [[COALESCENT], [COALESCENT, CONFIDENCE]].map((differences) =>
      variantSummary(input({ check: validCheck(differences) })),
    );

    expect(summaries).toStrictEqual([
      'Variant of "ebola/100, time tree". Changed: coalescent prior.',
      'Variant of "ebola/100, time tree". Changed: coalescent prior, confidence intervals.',
    ]);
  });

  test("another analysis is a change", () => {
    expect(variantSummary(input({ check: { ...validCheck([]), command: "clock" } }))).toBe(
      'Variant of "ebola/100, time tree". Changed: analysis.',
    );
  });

  test("on the page of its source the stale bar lists the changes, so the summary names only the source", () => {
    expect(
      [validCheck([]), validCheck([COALESCENT])].map((check) => variantSummary(input({ shownRunId: "source", check }))),
    ).toStrictEqual(["Variant of this run. No settings changed yet.", "Variant of this run"]);
  });

  test.each<[string, Partial<VariantInput>]>([
    ["the check is pending", { check: undefined }],
    ["the settings are invalid", { check: invalidCheck() }],
    ["the source is unknown", { source: undefined, check: validCheck([COALESCENT]) }],
    ["the page shows another run", { shownRunId: "other", check: validCheck([COALESCENT]) }],
  ])("without a comparison with the source (%s) the summary names only the source", (_case, changes) => {
    expect(variantSummary(input(changes))).toMatch(/^Variant of "[^"]+"$/u);
  });

  test("a draft that starts from no run and the read-only options of a run are no variant", () => {
    expect([
      variantSummary(input({ fromRunId: undefined, check: validCheck([COALESCENT]) })),
      variantSummary(input({ viewing: true, shownRunId: "source", check: validCheck([COALESCENT]) })),
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
    ...changes,
  };
}

function setting(key: string, label: string, purpose: SettingDifference["purpose"]): SettingDifference {
  return { kind: "setting", key, label, group: label, purpose, second: true };
}
