import type { IndexedConfig, IndexedDataset } from "@neherlab/app-contracts";

import {
  configInputs,
  exampleDataUrl,
  indexedDatasetInputs,
  localInputAssignment,
  uploadedConfigInputs,
} from "../exampleData";

const ZIKA_20: IndexedDataset = {
  name: "zika/20",
  inputs: [
    { kind: "tree", path: "zika/20/tree.nwk" },
    { kind: "alignment", path: "zika/20/aln.fasta.xz" },
    { kind: "metadata", path: "zika/20/metadata.tsv" },
  ],
};

const MUGRATION: IndexedConfig = {
  path: "zika/20/mugration.yaml",
  command: "mugration",
  title: "Mugration of the zika dataset",
  folder: "zika/20",
  content: "tree: tree.nwk\n",
  inputs: {
    metadata: ["zika/20/metadata.tsv"],
    tree: ["zika/20/tree.nwk"],
    weights: ["zika/country_weights.csv"],
  },
};

const TIMETREE: IndexedConfig = {
  path: "flu/h3n2/20/timetree.yaml",
  command: "timetree",
  title: "Timetree of the flu dataset",
  folder: "flu/h3n2/20",
  content: "tree: tree.nwk\n",
  inputs: { alignment: ["flu/h3n2/20/aln.fasta.xz"] },
};

describe("example data", () => {
  test.each([
    ["/", "index.json", "/data/index.json"],
    ["/treetime/", "index.json", "/treetime/data/index.json"],
    ["/treetime", "zika/20/tree.nwk", "/treetime/data/zika/20/tree.nwk"],
    ["./", "lassa/L/my tree#1.nwk", "./data/lassa/L/my%20tree%231.nwk"],
  ])("the base path %s puts %s at %s", (base, path, url) => {
    expect(exampleDataUrl(base, path)).toBe(url);
  });

  test("a dataset uploads the input files of the command, alignments as lists", () => {
    expect(indexedDatasetInputs(ZIKA_20, "timetree")).toStrictEqual([
      { key: "tree", list: false, path: "zika/20/tree.nwk" },
      { key: "metadata", list: false, path: "zika/20/metadata.tsv" },
      { key: "alignment", list: true, path: "zika/20/aln.fasta.xz" },
    ]);
  });

  test("a dataset uploads only the files of inputs the command has", () => {
    expect(indexedDatasetInputs(ZIKA_20, "mugration").map((input) => input.path)).toStrictEqual([
      "zika/20/tree.nwk",
      "zika/20/metadata.tsv",
    ]);
  });

  test.each([
    [
      "/home/user/.treetime/examples",
      [
        { key: "tree", value: "/home/user/.treetime/examples/zika/20/tree.nwk", label: "zika/20/tree.nwk" },
        { key: "metadata", value: "/home/user/.treetime/examples/zika/20/metadata.tsv", label: "zika/20/metadata.tsv" },
        {
          key: "alignment",
          value: ["/home/user/.treetime/examples/zika/20/aln.fasta.xz"],
          label: "zika/20/aln.fasta.xz",
        },
      ],
    ],
    [
      "C:\\TreeTime\\examples\\",
      [
        { key: "tree", value: "C:\\TreeTime\\examples\\zika\\20\\tree.nwk", label: "zika/20/tree.nwk" },
        { key: "metadata", value: "C:\\TreeTime\\examples\\zika\\20\\metadata.tsv", label: "zika/20/metadata.tsv" },
        { key: "alignment", value: ["C:\\TreeTime\\examples\\zika\\20\\aln.fasta.xz"], label: "zika/20/aln.fasta.xz" },
      ],
    ],
  ])("the desktop app fills a dataset from the examples folder %s", (folder, assignments) => {
    expect(indexedDatasetInputs(ZIKA_20, "timetree").map((input) => localInputAssignment(folder, input))).toStrictEqual(
      assignments,
    );
  });

  test("a config replaces each input it reads with the uploads, lists as lists", () => {
    const uploaded = [...configInputs(MUGRATION), ...configInputs(TIMETREE)].map((input, index) => ({
      input,
      stored: input.paths.map((path) => ({
        name: path.slice(path.lastIndexOf("/") + 1),
        path: `/up/${index}`,
        size: 7,
      })),
    }));

    expect(uploadedConfigInputs(uploaded, "r1")).toStrictEqual({
      replace_inputs: {
        metadata: "/up/0",
        tree: "/up/1",
        weights: "/up/2",
        alignment: ["/up/3"],
      },
      sources: {
        metadata: { label: "metadata.tsv", size: 7, origin: "upload", run_id: "r1" },
        tree: { label: "tree.nwk", size: 7, origin: "upload", run_id: "r1" },
        weights: { label: "country_weights.csv", size: 7, origin: "upload", run_id: "r1" },
        alignment: { label: "aln.fasta.xz", size: 7, origin: "upload", run_id: "r1" },
      },
    });
  });
});
