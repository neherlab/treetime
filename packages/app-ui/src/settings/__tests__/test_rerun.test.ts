import type { RunRecord } from "@neherlab/app-contracts";

import { commandSettings } from "../catalog";
import { defaultConfig } from "../config";
import { rerunDraft } from "../rerun";

const DEFAULT_CLOCK_CONFIG = defaultConfig(commandSettings("clock").settings);

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

describe("edit and run again", () => {
  test("keeps the changed settings and inputs of the run, drops its output paths and leaves its default outputs to the new settings", () => {
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
      selection: [],
      labels: { tree: "tree.nwk", metadata: "metadata.tsv" },
    });
  });

  test("keeps an output selection that the run chose", () => {
    const draft = rerunDraft({ ...RECORD, changed_settings: ["keep_root", "output_selection"] });

    expect(draft.config["output_selection"]).toStrictEqual(["nwk", "auspice"]);
  });

  test("leaves a setting that the run left at its default to the default of the new draft", () => {
    const draft = rerunDraft({ ...RECORD, config: { ...RECORD.config, allow_negative_rate: true } });

    expect(draft.config["allow_negative_rate"]).toStrictEqual(DEFAULT_CLOCK_CONFIG["allow_negative_rate"]);
  });

  test("the original record is not changed", () => {
    const before = JSON.stringify(RECORD);
    rerunDraft(RECORD);

    expect(JSON.stringify(RECORD)).toStrictEqual(before);
  });
});
