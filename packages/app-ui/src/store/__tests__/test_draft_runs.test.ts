import type { UiDraft } from "@neherlab/app-contracts";

import { draftInputPaths, forgetRuns, freshDraft, referencedRuns } from "../draftSchema";

const FRESH = freshDraft("timetree");

const DRAFT: UiDraft = {
  ...FRESH,
  config: {
    ...FRESH.config,
    tree: "/treetime/runs/r1/inputs/tree.nwk",
    alignment: ["/treetime/runs/r1/inputs/aln.fasta"],
    metadata: "/treetime/runs/r2/inputs/dates.tsv",
    topology_order_list: "/data/order.txt",
    clock_rate: 0.001,
  },
  sources: {
    tree: { label: "tree.nwk", origin: { kind: "upload" }, size: 10, run_id: "r1" },
    alignment: { label: "aln.fasta", origin: { kind: "upload" }, size: 20, run_id: "r1" },
    metadata: { label: "dates.tsv", origin: { kind: "run" }, run_id: "r2" },
    topology_order_list: { label: "order.txt", origin: { kind: "local" } },
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
        topology_order_list: "/data/order.txt",
        clock_rate: 0.001,
      },
      sources: {
        metadata: { label: "dates.tsv", origin: { kind: "run" }, run_id: "r2" },
        topology_order_list: { label: "order.txt", origin: { kind: "local" } },
      },
      from_run_id: "r3",
      search: "rate",
    };

    const lost = {
      runs: new Set(["r1"]),
      paths: new Set(["/treetime/runs/r1/inputs/tree.nwk", "/treetime/runs/r1/inputs/aln.fasta"]),
    };

    expect(forgetRuns(DRAFT, lost)).toStrictEqual({
      draft: expected,
      dropped: ["tree.nwk", "aln.fasta"],
      examples: [],
    });
  });

  test("forgetting the earlier run drops the form's link to it and no file", () => {
    const { from_run_id: _fromRun, ...expected } = DRAFT;

    expect(forgetRuns(DRAFT, { runs: new Set(["r3"]), paths: new Set() })).toStrictEqual({
      draft: expected,
      dropped: [],
      examples: [],
    });
  });

  test("forgetting no run keeps the draft", () => {
    expect(forgetRuns(DRAFT, { runs: new Set(), paths: new Set() })).toStrictEqual({
      draft: DRAFT,
      dropped: [],
      examples: [],
    });
  });

  test("an input of a deleted run whose stored copy another run holds stays in the draft", () => {
    const { from_run_id: _fromRun, ...expected } = DRAFT;

    expect(forgetRuns(DRAFT, { runs: new Set(["r2", "r3"]), paths: new Set() })).toStrictEqual({
      draft: expected,
      dropped: [],
      examples: [],
    });
  });

  test("a stored input that is gone drops its file also when its run exists", () => {
    const forgotten = forgetRuns(DRAFT, { runs: new Set(), paths: new Set(["/treetime/runs/r2/inputs/dates.tsv"]) });

    expect([Object.keys(forgotten.draft.sources), forgotten.dropped]).toStrictEqual([
      ["tree", "alignment", "topology_order_list"],
      ["dates.tsv"],
    ]);
  });

  test("the input paths of a draft are the paths of its sources, once each", () => {
    expect(draftInputPaths(DRAFT)).toStrictEqual([
      "/treetime/runs/r1/inputs/tree.nwk",
      "/treetime/runs/r1/inputs/aln.fasta",
      "/treetime/runs/r2/inputs/dates.tsv",
      "/data/order.txt",
    ]);
  });
});
