import { describe, expect, test } from "vitest";

import { zPickedFiles, zPickFilesRequest } from "../files";

describe("file picker request", () => {
  test("accepts a complete request", () => {
    const request = { title: "Tree", extensions: ["nwk", "nexus"], multiple: true };

    expect(zPickFilesRequest.parse(request)).toStrictEqual(request);
  });

  test.each([
    { name: "an unknown field", value: { title: "Tree", extensions: [], multiple: false, path: "/etc" } },
    { name: "a missing field", value: { title: "Tree", extensions: [] } },
    { name: "a non-string extension", value: { title: "Tree", extensions: [1], multiple: false } },
    { name: "a non-boolean multiple", value: { title: "Tree", extensions: [], multiple: "yes" } },
  ])("rejects $name", ({ value }) => {
    expect(zPickFilesRequest.safeParse(value).success).toBe(false);
  });
});

describe("picked files", () => {
  test("accepts a list of paths", () => {
    expect(zPickedFiles.parse(["/data/tree.nwk"])).toStrictEqual(["/data/tree.nwk"]);
  });

  test("rejects a non-string path", () => {
    expect(zPickedFiles.safeParse(["/data/tree.nwk", 3]).success).toBe(false);
  });
});
