import type { UiDraft } from "@neherlab/app-contracts";

import { forgetRuns, freshDraft, referencedRuns } from "../draftSchema";

const FRESH = freshDraft("timetree");

const DRAFT: UiDraft = {
  ...FRESH,
  config: {
    ...FRESH.config,
    tree: "/treetime/runs/r1/inputs/tree.nwk",
    alignment: ["/treetime/runs/r1/inputs/aln.fasta"],
    metadata: "/treetime/runs/r2/inputs/dates.tsv",
    custom_gtr: "/data/gtr.json",
    clock_rate: 0.001,
  },
  sources: {
    tree: { label: "tree.nwk", origin: "upload", size: 10, run_id: "r1" },
    alignment: { label: "aln.fasta", origin: "upload", size: 20, run_id: "r1" },
    metadata: { label: "dates.tsv", origin: "run", run_id: "r2" },
    custom_gtr: { label: "gtr.json", origin: "local" },
  },
  upload_run_id: "r1",
  from_run_id: "r3",
  search: "rate",
};

describe("runs of a draft", () => {
  test("the referenced runs are the upload run, the earlier run, and the runs of the sources, once each", () => {
    expect(referencedRuns(DRAFT)).toStrictEqual(["r1", "r3", "r2"]);
  });

  test("a fresh draft references no run", () => {
    expect(referencedRuns(FRESH)).toStrictEqual([]);
  });

  test("forgetting the upload run resets its files to their defaults and keeps every other setting", () => {
    const expected: UiDraft = {
      ...FRESH,
      config: {
        ...FRESH.config,
        metadata: "/treetime/runs/r2/inputs/dates.tsv",
        custom_gtr: "/data/gtr.json",
        clock_rate: 0.001,
      },
      sources: {
        metadata: { label: "dates.tsv", origin: "run", run_id: "r2" },
        custom_gtr: { label: "gtr.json", origin: "local" },
      },
      from_run_id: "r3",
      search: "rate",
    };

    expect(forgetRuns(DRAFT, new Set(["r1"]))).toStrictEqual({ draft: expected, dropped: ["tree.nwk", "aln.fasta"] });
  });

  test("forgetting the earlier run drops the form's link to it and no file", () => {
    const { from_run_id: _fromRun, ...expected } = DRAFT;

    expect(forgetRuns(DRAFT, new Set(["r3"]))).toStrictEqual({ draft: expected, dropped: [] });
  });

  test("forgetting no run keeps the draft", () => {
    expect(forgetRuns(DRAFT, new Set())).toStrictEqual({ draft: DRAFT, dropped: [] });
  });
});
