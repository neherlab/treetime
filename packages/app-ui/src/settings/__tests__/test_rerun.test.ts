import type { RunRecord } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { rerunDraft } from "../rerun";

const RECORD: RunRecord = {
  id: "r1",
  title: "Baseline",
  command: "clock",
  config: {
    tree: "/runs/r0/inputs/tree.nwk",
    metadata: "/data/zika/86/metadata.tsv",
    keep_root: true,
    output_all: "/runs/r1/out",
    output_selection: ["Nwk", "Auspice"],
  },
  status: "ok",
  pinned: false,
  created_at: "2026-09-25T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [
    { setting: "tree", path: "/runs/r0/inputs/tree.nwk", size: 10, sha256: "a" },
    { setting: "metadata", path: "/data/zika/86/metadata.tsv", size: 20, sha256: "b" },
  ],
  changed_settings: ["keep_root"],
  headline: {},
  output_files: [],
};

describe("edit and run again", () => {
  test("keeps the settings and inputs of the run and drops its output paths", () => {
    const draft = rerunDraft(RECORD);

    expect({
      tree: draft.config["tree"],
      keepRoot: draft.config["keep_root"],
      outputAll: draft.config["output_all"],
      selection: draft.config["output_selection"],
      labels: draft.inputLabels,
      title: draft.title,
    }).toStrictEqual({
      tree: "/runs/r0/inputs/tree.nwk",
      keepRoot: true,
      outputAll: null,
      selection: ["Nwk", "Auspice"],
      labels: { tree: "tree.nwk", metadata: "metadata.tsv" },
      title: "Baseline (edited)",
    });
  });

  test("the original record is not changed", () => {
    const before = JSON.stringify(RECORD);
    rerunDraft(RECORD);

    expect(JSON.stringify(RECORD)).toStrictEqual(before);
  });
});
