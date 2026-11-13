import type { DetectedFile, InputSlot } from "@neherlab/app-contracts";

import { placeFiles } from "../filePlacement";

const TREE = slot("tree", "Tree", false);

const SLOTS: InputSlot[] = [TREE, slot("alignment", "Alignment", true), slot("metadata", "Metadata", false)];

function slot(kind: InputSlot["kind"], label: string, list: boolean): InputSlot {
  return { kind, label, list, need: "required", formats: "", extensions: [] };
}

function file(path: string, kind: DetectedFile["kind"], reason?: string): DetectedFile {
  const detected: DetectedFile = { path, kind, commands: [] };

  if (reason !== undefined) {
    detected.reason = reason;
  }

  return detected;
}

describe("placement of dropped files", () => {
  test("each data file fills the input of its kind, alignments as a list", () => {
    const detected = [
      file("/up/aln1.fasta", "alignment"),
      file("/up/tree.nwk", "tree"),
      file("/up/meta.tsv", "metadata"),
      file("/up/aln2.fasta", "alignment"),
    ];

    expect(placeFiles(detected, SLOTS, "Time tree")).toStrictEqual({
      inputs: [
        { key: "tree", list: false, files: [1] },
        { key: "alignment", list: true, files: [0, 3] },
        { key: "metadata", list: false, files: [2] },
      ],
      configs: [],
      packages: [],
      unplaced: [],
    });
  });

  test("configurations and run packages are kept apart, and unknown and surplus files are listed with the reason", () => {
    const detected = [
      file("/up/a.nwk", "tree"),
      file("/up/b.nwk", "tree"),
      file("/up/run.zip", "run-package"),
      file("/up/notes.txt", "unknown", "not a tree: no opening parenthesis"),
      file("/up/one.yaml", "configuration"),
      file("/up/two.yaml", "configuration"),
    ];

    expect(placeFiles(detected, SLOTS, "Time tree")).toStrictEqual({
      inputs: [{ key: "tree", list: false, files: [0] }],
      configs: [4],
      packages: [2],
      unplaced: [
        { index: 1, reason: "The analysis reads one tree; a.nwk fills it" },
        { index: 3, reason: "not a tree: no opening parenthesis" },
        { index: 5, reason: "Only one configuration loads at a time; one.yaml loads" },
      ],
    });
  });

  test("a data file whose kind the analysis does not read is listed with the reason", () => {
    const detected = [file("/up/aln.fasta", "alignment")];

    expect(placeFiles(detected, [TREE], "Clock signal").unplaced).toStrictEqual([
      { index: 0, reason: "The Clock signal analysis reads no alignment" },
    ]);
  });
});
