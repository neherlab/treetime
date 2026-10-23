import { commandSettings } from "../catalog";
import { defaultConfig } from "../config";
import { folderName, inputFactsRequest, pathList, runInputAssignments, slotFactsText, slotProblems } from "../inputs";
import { withSetting } from "../json";

describe("inputs", () => {
  test("the inputs of an earlier run become settings with their recorded paths", () => {
    const inputs = [
      { setting: "tree", path: "/runs/a/inputs/tree.nwk", size: 1, sha256: "x" },
      { setting: "alignment", path: "/runs/a/inputs/one.fasta", size: 1, sha256: "y" },
      { setting: "alignment", path: "/runs/a/inputs/two.fasta", size: 1, sha256: "z" },
    ];

    expect(runInputAssignments("ancestral", inputs)).toStrictEqual([
      { key: "tree", value: "/runs/a/inputs/tree.nwk", label: "tree.nwk" },
      {
        key: "alignment",
        value: ["/runs/a/inputs/one.fasta", "/runs/a/inputs/two.fasta"],
        label: "one.fasta, two.fasta",
      },
    ]);
  });

  test("the input check sends the command and its configuration", () => {
    const config = withSetting(defaultConfig(commandSettings("mugration").settings), "tree", "t.nwk");

    expect(inputFactsRequest("mugration", config)).toStrictEqual({ command: "mugration", config });
  });

  test("a path setting lists its non-empty paths", () => {
    expect([pathList("t.nwk"), pathList(["a.fasta", ""]), pathList(undefined), pathList("")]).toStrictEqual([
      ["t.nwk"],
      ["a.fasta"],
      [],
      [],
    ]);
  });

  test("no input means no input check", () => {
    expect(inputFactsRequest("prune", defaultConfig(commandSettings("prune").settings))).toStrictEqual(null);
  });
});

describe("slot facts", () => {
  const facts = {
    tree: { tips: 86, internal_nodes: 61, polytomies: 15, unnamed_tips: 0, duplicate_node_names: [] },
    alignments: [
      {
        number: 1,
        path: "zika/aln.fasta",
        sequences: 86,
        min_length: 10807,
        max_length: 10807,
        duplicate_names: [],
      },
    ],
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
      ["86 tips, 61 internal nodes, 15 polytomies"],
      ["86 sequences, 10,807 sites"],
      ["86 rows; ID column name; date column date"],
      ["86 rows; ID column name"],
    ]);
  });

  test("an unequal alignment says it is not aligned", () => {
    const uneven = {
      ...facts,
      alignments: [{ number: 1, path: "a.fasta", sequences: 2, min_length: 9, max_length: 10, duplicate_names: [] }],
    };

    expect(slotFactsText("alignment", uneven, false)).toStrictEqual([
      "2 sequences, lengths 9 to 10, so it is not aligned",
    ]);
  });

  test("several alignment files show their facts and problems per file", () => {
    const several = {
      ...facts,
      alignments: [
        {
          number: 1,
          path: "/runs/r1/inputs/1/aligned.fasta",
          sequences: 86,
          min_length: 9,
          max_length: 9,
          duplicate_names: [],
        },
        {
          number: 2,
          path: "/runs/r1/inputs/2/aligned.fasta",
          sequences: 80,
          min_length: 4,
          max_length: 5,
          duplicate_names: [],
        },
      ],
      problems: [{ input: "alignment" as const, alignment: 3, message: "When reading file 'mp.fasta': no such file" }],
    };

    expect([slotFactsText("alignment", several, false), slotProblems("alignment", several)]).toStrictEqual([
      [
        "Alignment 1 (aligned.fasta): 86 sequences, 9 sites",
        "Alignment 2 (aligned.fasta): 80 sequences, lengths 4 to 5, so it is not aligned",
      ],
      ["Alignment 3: When reading file 'mp.fasta': no such file"],
    ]);
  });

  test("a slot shows the problems of its own input", () => {
    expect([slotProblems("tree", facts), slotProblems("metadata", facts)]).toStrictEqual([["bad"], []]);
  });
});

describe("folder name of a path", () => {
  test.each([
    ["/home/user/zika/ancestral.yaml", "/home/user/zika"],
    ["C:\\data\\zika\\ancestral.yaml", "C:\\data\\zika"],
    ["/ancestral.yaml", "/"],
    ["ancestral.yaml", undefined],
  ])("%s is in %s", (path, expected) => {
    expect(folderName(path)).toBe(expected);
  });
});
