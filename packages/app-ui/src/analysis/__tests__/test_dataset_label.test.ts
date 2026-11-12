import type { UiDraft } from "@neherlab/app-contracts";

import { freshDraft } from "../../store/draftSchema";
import { datasetLabel } from "../datasetLabel";

const FRESH = freshDraft("timetree");

describe("dataset label of the workspace header", () => {
  test("an example dataset shows its name", () => {
    const draft: UiDraft = {
      ...FRESH,
      sources: {
        metadata: { label: "metadata.tsv", origin: { kind: "example", path: "ebola/100/metadata.tsv" }, run_id: "u1" },
        tree: { label: "tree.nwk", origin: { kind: "example", path: "ebola/100/tree.nwk" }, run_id: "u1" },
      },
    };

    expect(datasetLabel(draft)).toBe("ebola/100");
  });

  test("an example dataset of the data folder on the computer shows its name", () => {
    const draft: UiDraft = {
      ...FRESH,
      sources: { tree: { label: "zika/86/tree.nwk", origin: { kind: "dataset" } } },
    };

    expect(datasetLabel(draft)).toBe("zika/86");
  });

  test("own data shows the file names in the order of the input slots", () => {
    const draft: UiDraft = {
      ...FRESH,
      sources: {
        metadata: { label: "dates.tsv", origin: { kind: "upload" }, run_id: "u1" },
        tree: { label: "my.nwk", origin: { kind: "upload" }, run_id: "u1" },
      },
    };

    expect(datasetLabel(draft)).toBe("my.nwk, dates.tsv");
  });

  test("an example tree with own metadata shows the file names", () => {
    const draft: UiDraft = {
      ...FRESH,
      sources: {
        tree: { label: "tree.nwk", origin: { kind: "example", path: "ebola/100/tree.nwk" }, run_id: "u1" },
        metadata: { label: "dates.tsv", origin: { kind: "upload" }, run_id: "u1" },
      },
    };

    expect(datasetLabel(draft)).toBe("tree.nwk, dates.tsv");
  });

  test("a draft without data has no label", () => {
    expect(datasetLabel(FRESH)).toBeUndefined();
  });
});
