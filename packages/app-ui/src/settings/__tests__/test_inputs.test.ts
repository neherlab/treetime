import { describe, expect, test } from "vitest";

import { defaultConfig } from "../config";
import {
  datasetInputs,
  filledSlots,
  inputFactsRequest,
  pathList,
  runInputAssignments,
  slotFactsText,
  slotProblem,
} from "../inputs";
import { setAt } from "../json";
import { commandSettings } from "../schema";

const ZIKA_86 = { name: "zika/86", files: ["aln.fasta.xz", "metadata.tsv", "tree.nwk", "zika.phylip.xz"] };

describe("inputs", () => {
  test("a dataset fills the slots of the command from the data directory", () => {
    expect(datasetInputs("data", ZIKA_86, "timetree")).toStrictEqual([
      { key: "tree", value: "data/zika/86/tree.nwk", label: "zika/86/tree.nwk" },
      { key: "metadata", value: "data/zika/86/metadata.tsv", label: "zika/86/metadata.tsv" },
      { key: "alignment", value: ["data/zika/86/aln.fasta.xz"], label: "zika/86/aln.fasta.xz" },
    ]);
  });

  test("a dataset fills only the slots the command has", () => {
    expect(datasetInputs("/srv/data", ZIKA_86, "mugration").map((input) => input.value)).toStrictEqual([
      "/srv/data/zika/86/tree.nwk",
      "/srv/data/zika/86/metadata.tsv",
    ]);
  });

  test("the inputs of an earlier run become settings with their recorded paths", () => {
    const inputs = [
      { setting: "tree", path: "/runs/a/inputs/tree.nwk", size: 1, sha256: "x" },
      { setting: "alignment", path: "/runs/a/inputs/one.fasta", size: 1, sha256: "y" },
      { setting: "alignment", path: "/runs/a/inputs/two.fasta", size: 1, sha256: "z" },
    ];

    expect(runInputAssignments(inputs)).toStrictEqual([
      { key: "tree", value: "/runs/a/inputs/tree.nwk", label: "tree.nwk" },
      {
        key: "alignment",
        value: ["/runs/a/inputs/one.fasta", "/runs/a/inputs/two.fasta"],
        label: "one.fasta, two.fasta",
      },
    ]);
  });

  test("filled slots are the inputs with a path", () => {
    const config = setAt(defaultConfig(commandSettings("clock").specs), ["tree"], "t.nwk");

    expect([...filledSlots("clock", config)]).toStrictEqual(["tree"]);
  });

  test("the input check reads the inputs and metadata settings of the command", () => {
    const config = setAt(
      setAt(defaultConfig(commandSettings("mugration").specs), ["tree"], "t.nwk"),
      ["metadata"],
      "m.tsv",
    );

    expect(inputFactsRequest("mugration", config)).toStrictEqual({
      tree: "t.nwk",
      metadata: "m.tsv",
      alignment: [],
      metadata_id_columns: ["strain", "name", "accession"],
      metadata_delimiters: [",", "\t", ";"],
      date_column: null,
    });
  });

  test("a path setting lists its non-empty paths", () => {
    expect([pathList("t.nwk"), pathList(["a.fasta", ""]), pathList(null), pathList("")]).toStrictEqual([
      ["t.nwk"],
      ["a.fasta"],
      [],
      [],
    ]);
  });

  test("no input means no input check", () => {
    expect(inputFactsRequest("prune", defaultConfig(commandSettings("prune").specs))).toStrictEqual(null);
  });
});

describe("slot facts", () => {
  const facts = {
    tree: { tips: 86, internal_nodes: 61, polytomies: 15, unnamed_tips: 0, duplicate_tip_names: [] },
    alignment: { sequences: 86, min_length: 10807, max_length: 10807, duplicate_names: [] },
    metadata: { rows: 86, columns: ["name", "date"], id_column: "name", date_column: "date" },
    problems: [{ input: "tree" as const, message: "bad" }],
  };

  test("each slot summarizes its facts", () => {
    expect([
      slotFactsText("tree", facts, true),
      slotFactsText("alignment", facts, true),
      slotFactsText("metadata", facts, true),
      slotFactsText("metadata", facts, false),
    ]).toStrictEqual([
      "86 tips, 61 internal nodes, 15 polytomies",
      "86 sequences, 10,807 sites",
      "86 rows; ID column name; date column date",
      "86 rows; ID column name",
    ]);
  });

  test("an unequal alignment says it is not aligned", () => {
    const uneven = { ...facts, alignment: { sequences: 2, min_length: 9, max_length: 10, duplicate_names: [] } };

    expect(slotFactsText("alignment", uneven, false)).toStrictEqual(
      "2 sequences, lengths 9 to 10, so it is not aligned",
    );
  });

  test("a slot shows the problem of its own input", () => {
    expect([slotProblem("tree", facts), slotProblem("metadata", facts)]).toStrictEqual(["bad", null]);
  });
});
