import { describe, expect, test } from "vitest";

import { dropzoneAccept } from "../analysis/fileAccept";

describe("file dialog accept filter", () => {
  test("extensions become dotted patterns under the application wildcard type", () => {
    expect(dropzoneAccept(["nwk", "fasta.xz"])).toStrictEqual({ accept: { "application/*": [".nwk", ".fasta.xz"] } });
  });

  test("a setting without extensions sets no accept filter", () => {
    expect(dropzoneAccept([])).toStrictEqual({});
  });
});
