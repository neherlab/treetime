import type { UploadAllowance } from "@neherlab/app-contracts";

import { oversizedUpload, storeWithinAllowance } from "../uploadAllowance";

const ALLOWANCE: UploadAllowance = { remaining: 1000, limit: 1024, limit_text: "1 KiB" };

describe("upload allowance", () => {
  test("files that fit the remaining allowance together pass", () => {
    expect(
      oversizedUpload(
        [
          { name: "tree.nwk", size: 600 },
          { name: "metadata.tsv", size: 400 },
        ],
        ALLOWANCE,
      ),
    ).toBeUndefined();
  });

  test("the first file that takes the total over the allowance is named with the limit", () => {
    expect(
      oversizedUpload(
        [
          { name: "tree.nwk", size: 600 },
          { name: "aln.fasta", size: 401 },
          { name: "metadata.tsv", size: 1 },
        ],
        ALLOWANCE,
      ),
    ).toBe("the inputs of a run are limited to 1 KiB in total; `aln.fasta` does not fit");
  });

  test("a file over the allowance is refused before any transfer", async () => {
    const transfers: string[][] = [];

    const store = (files: readonly { name: string; size: number }[]) => {
      transfers.push(files.map(({ name }) => name));

      return Promise.resolve([]);
    };

    await expect(
      storeWithinAllowance([{ name: "aln.fasta", size: 2048 }], () => Promise.resolve(ALLOWANCE), store),
    ).rejects.toThrow("the inputs of a run are limited to 1 KiB in total; `aln.fasta` does not fit");
    expect(transfers).toStrictEqual([]);
  });

  test("files within the allowance are stored", async () => {
    const stored = await storeWithinAllowance(
      [{ name: "tree.nwk", size: 10 }],
      () => Promise.resolve(ALLOWANCE),
      (files) => Promise.resolve(files.map(({ name, size }) => ({ name, path: `/inputs/${name}`, size }))),
    );

    expect(stored).toStrictEqual([{ name: "tree.nwk", path: "/inputs/tree.nwk", size: 10 }]);
  });
});
