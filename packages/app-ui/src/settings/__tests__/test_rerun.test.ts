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
    output_selection: ["nwk", "auspice"],
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
  warnings: [],
};

const ANCESTRAL: RunRecord = {
  id: "r2",
  title: "Ancestral zika",
  command: "ancestral",
  config: {
    tree: "/data/zika/86/tree.nwk",
    alignment: ["/data/zika/86/aln.fasta.xz"],
    model: "jc69",
    gap_fill: "all",
    reconstruct_tip_states: true,
    output_all: "/runs/r2/out",
    output_selection: ["nwk", "augur-node-data", "gtr", "auspice"],
  },
  status: "ok",
  pinned: false,
  created_at: "2026-09-25T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [
    { setting: "tree", path: "/data/zika/86/tree.nwk", size: 10, sha256: "a" },
    { setting: "alignment", path: "/data/zika/86/aln.fasta.xz", size: 20, sha256: "b" },
  ],
  changed_settings: ["model", "gap_fill", "reconstruct_tip_states"],
  headline: {},
  output_files: [],
  warnings: [],
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
    }).toStrictEqual({
      tree: "/runs/r0/inputs/tree.nwk",
      keepRoot: true,
      outputAll: undefined,
      selection: ["nwk", "auspice"],
      labels: { tree: "tree.nwk", metadata: "metadata.tsv" },
    });
  });

  test("the original record is not changed", () => {
    const before = JSON.stringify(RECORD);
    rerunDraft(RECORD);

    expect(JSON.stringify(RECORD)).toStrictEqual(before);
  });
});

describe("analyze homoplasy from an ancestral run", () => {
  test("carries the inputs and the changed settings that homoplasy also has", () => {
    const draft = rerunDraft(ANCESTRAL, "homoplasy");

    expect({
      tree: draft.config["tree"],
      alignment: draft.config["alignment"],
      model: draft.config["model"],
      gapFill: draft.config["gap_fill"],
      tipStates: draft.config["reconstruct_tip_states"],
      outputAll: draft.config["output_all"],
      selection: draft.config["output_selection"],
      labels: draft.inputLabels,
    }).toStrictEqual({
      tree: "/data/zika/86/tree.nwk",
      alignment: ["/data/zika/86/aln.fasta.xz"],
      model: "jc69",
      gapFill: "all",
      tipStates: undefined,
      outputAll: undefined,
      selection: [],
      labels: { tree: "tree.nwk", alignment: "aln.fasta.xz" },
    });
  });
});
