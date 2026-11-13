import type { DetectedFile } from "@neherlab/app-contracts";

import { detectedFiles, fileHead, withRoles } from "../openedFiles";

describe("opened files", () => {
  test("the head of a small file holds the whole file", async () => {
    expect(await fileHead(new Blob(["tree: tree.nwk\n"]))).toStrictEqual({ head: "tree: tree.nwk\n", complete: true });
  });

  test("the head of a large file holds its first bytes only", async () => {
    const head = await fileHead(new Blob([">A\n", "A".repeat(100_000)]));

    expect([head.head.length, head.complete]).toStrictEqual([65_536, false]);
  });

  test("only the inputs among the opened files are read as inputs, in the order of the files", () => {
    const files = withRoles(
      [{ name: "run.zip" }, { name: "tree.nwk" }, { name: "run.yaml" }, { name: "aln.fasta" }],
      ["run-package", "input", "configuration", "input"],
    );

    const inputs: DetectedFile[] = [
      { path: "/runs/r1/inputs/1/tree.nwk", kind: "tree", commands: [] },
      { path: "/runs/r1/inputs/2/aln.fasta", kind: "alignment", commands: [] },
    ];

    expect([
      files.filter(({ role }) => role === "input").map(({ index }) => index),
      detectedFiles(files, inputs).map(({ path, kind }) => [path, kind]),
    ]).toStrictEqual([
      [1, 3],
      [
        ["run.zip", "run-package"],
        ["/runs/r1/inputs/1/tree.nwk", "tree"],
        ["run.yaml", "configuration"],
        ["/runs/r1/inputs/2/aln.fasta", "alignment"],
      ],
    ]);
  });

  test("a package and a configuration alone leave no input to read", () => {
    const files = withRoles([{ name: "run.zip" }, { name: "run.yaml" }], ["run-package", "configuration"]);

    expect(files.filter(({ role }) => role === "input")).toStrictEqual([]);
  });
});
