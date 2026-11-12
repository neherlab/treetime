import type { UiDraftSource } from "@neherlab/app-contracts";

import { draftExample } from "../draftExample";

const EBOLA_TREE: UiDraftSource = {
  label: "tree.nwk",
  origin: { kind: "example", path: "ebola/100/tree.nwk" },
  run_id: "u1",
};

const EBOLA_METADATA: UiDraftSource = {
  label: "metadata.tsv",
  origin: { kind: "example", path: "ebola/100/metadata.tsv" },
  run_id: "u1",
};

const ZIKA_METADATA: UiDraftSource = {
  label: "metadata.tsv",
  origin: { kind: "example", path: "zika/20/metadata.tsv" },
  run_id: "u1",
};

const OWN_METADATA: UiDraftSource = { label: "dates.tsv", origin: { kind: "upload" }, run_id: "u1" };

const RUN_TREE: UiDraftSource = { label: "tree.nwk", origin: { kind: "run" }, run_id: "r1" };

describe("example dataset of a draft", () => {
  test.each<[string, UiDraftSource[], string | undefined, string | undefined]>([
    ["the files of one example name it", [EBOLA_TREE, EBOLA_METADATA], undefined, "ebola/100"],
    [
      "a file of the data folder on the computer names its folder",
      [{ label: "zika/86/tree.nwk", origin: { kind: "dataset" } }],
      undefined,
      "zika/86",
    ],
    ["own files name no example", [OWN_METADATA], "ebola/100", undefined],
    ["an example tree with own metadata names no example", [EBOLA_TREE, OWN_METADATA], "ebola/100", undefined],
    ["files of two examples name no example", [EBOLA_TREE, ZIKA_METADATA], undefined, undefined],
    ["the inputs of a run keep the example of the run", [RUN_TREE], "ebola/100", "ebola/100"],
    ["the inputs of a run of own files name no example", [RUN_TREE], undefined, undefined],
    ["the inputs of a run with an own file name no example", [RUN_TREE, OWN_METADATA], "ebola/100", undefined],
    ["a draft without files names no example", [], "ebola/100", undefined],
  ])("%s", (_name, sources, previous, expected) => {
    expect(draftExample(sources, previous)).toBe(expected);
  });
});
