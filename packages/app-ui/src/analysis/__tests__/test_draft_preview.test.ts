import type { InputFacts, RegressionStatistic } from "@neherlab/app-contracts";

import { inputFactRows, statisticText } from "../draftPreview";

const FACTS: InputFacts = {
  tree: { tips: 20, internal_nodes: 17, polytomies: 1, unnamed_tips: 0, duplicate_node_names: [], tip_names: [] },
  alignments: [
    {
      number: 1,
      path: "/runs/.inputs/a/aln.fasta",
      sequences: 20,
      min_length: 10807,
      max_length: 10807,
      duplicate_names: [],
    },
  ],
  metadata: {
    rows: 21,
    columns: ["strain", "date"],
    id_column: "strain",
    date_column: "date",
    dates: {
      readable: 19,
      unreadable: ["a"],
      exact_days: 19,
      on_day_1_or_15: 0,
      span: { earliest: { year: 2013.5, date: "2013-07-02" }, latest: { year: 2016.1, date: "2016-02-07" } },
    },
  },
  tips_without_metadata: ["t1"],
  detected: {},
  problems: [],
};

describe("input facts of a draft", () => {
  test("lists the tree, the alignment, the metadata, and the dates", () => {
    expect(inputFactRows(FACTS)).toStrictEqual([
      { label: "Tree", value: "20 tips, 17 internal nodes, 1 polytomy" },
      { label: "Alignment", value: "aln.fasta: 20 sequences, length 10807" },
      { label: "Metadata", value: "21 rows, ID column strain" },
      { label: "Dates", value: "19 readable, 1 unreadable, 2013-07-02 to 2016-02-07" },
      { label: "Tips without metadata", value: "1" },
    ]);
  });

  test("a draft with only a tree lists the tree", () => {
    const { metadata: _metadata, tips_without_metadata: _tips, ...treeOnly } = FACTS;

    expect(inputFactRows({ ...treeOnly, alignments: [] })).toStrictEqual([
      { label: "Tree", value: "20 tips, 17 internal nodes, 1 polytomy" },
    ]);
  });
});

describe("regression statistic text", () => {
  test.each<[RegressionStatistic, string]>([
    [
      {
        statistic: "clock-rate",
        label: "Clock rate",
        definition: "",
        unit: "/site/yr",
        value: { kind: "number", value: 0.000_765_4 },
      },
      "7.65e-4 /site/yr",
    ],
    [{ statistic: "r-squared", label: "R²", definition: "", value: { kind: "number", value: 0.912_34 } }, "0.912"],
    [
      {
        statistic: "root-date",
        label: "Root date",
        definition: "",
        value: { kind: "date", date: { year: 2013.2, date: "2013-03-14" } },
      },
      "2013-03-14",
    ],
  ])("formats %#", (statistic, expected) => {
    expect(statisticText(statistic)).toBe(expected);
  });
});
