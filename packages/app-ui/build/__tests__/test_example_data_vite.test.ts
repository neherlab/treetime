import type { DataIndex } from "@neherlab/app-contracts";

import { exampleDataTargets, parseDataIndex } from "../example-data-vite";

const INDEX: DataIndex = {
  files: [
    { path: "README.md", size: 10 },
    { path: "zika/20/tree.nwk", size: 20 },
    { path: "zika/20/metadata.tsv", size: 30 },
  ],
  datasets: [
    {
      name: "zika/20",
      inputs: [
        { kind: "tree", path: "zika/20/tree.nwk" },
        { kind: "metadata", path: "zika/20/metadata.tsv" },
      ],
    },
  ],
  configs: [],
  pipelines: [],
};

describe("example data plugin", () => {
  test("copies and serves the index and every listed file under data/, and nothing else", () => {
    expect(exampleDataTargets("/repo/data", INDEX)).toStrictEqual([
      { src: "/repo/data/index.json", dest: "data", rename: { stripBase: true } },
      { src: "/repo/data/README.md", dest: "data", rename: { stripBase: true } },
      { src: "/repo/data/zika/20/tree.nwk", dest: "data/zika/20", rename: { stripBase: true } },
      { src: "/repo/data/zika/20/metadata.tsv", dest: "data/zika/20", rename: { stripBase: true } },
    ]);
  });

  test("escapes glob characters in file names, so each target names exactly one file", () => {
    const index = { ...INDEX, files: [{ path: "flu/[draft] tree (1).nwk", size: 1 }] };

    expect(exampleDataTargets("/repo/data", index).map((target) => target.src)).toStrictEqual([
      "/repo/data/index.json",
      "/repo/data/flu/\\[draft\\] tree \\(1\\).nwk",
    ]);
  });

  test("reads an index whose files exist", () => {
    expect(parseDataIndex(JSON.stringify(INDEX), () => true)).toStrictEqual(INDEX);
  });

  test("rejects an index that does not match the schema", () => {
    const withoutPipelines = Object.fromEntries(Object.entries(INDEX).filter(([key]) => key !== "pipelines"));

    expect(() => parseDataIndex(JSON.stringify(withoutPipelines), () => true)).toThrow(
      "The example data index does not match the DataIndex schema",
    );
  });

  test("rejects an index that lists a missing file", () => {
    expect(() => parseDataIndex(JSON.stringify(INDEX), (path) => path !== "zika/20/metadata.tsv")).toThrow(
      "The example data index lists missing files: zika/20/metadata.tsv. Regenerate it with: just data-reindex",
    );
  });

  test.each(["../secret.txt", "/etc/passwd", "zika/../../secret.txt", "zika/./tree.nwk", ""])(
    "rejects an index that lists the path %p outside the data folder",
    (path) => {
      const index = { ...INDEX, files: [{ path, size: 1 }] };

      expect(() => parseDataIndex(JSON.stringify(index), () => true)).toThrow(
        `The example data index lists paths outside the data folder: ${path}`,
      );
    },
  );
});
