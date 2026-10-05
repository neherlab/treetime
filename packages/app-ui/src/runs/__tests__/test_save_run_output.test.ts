import { createApiClient } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";

import { saveRunOutput, type SaveTarget } from "../saveRunOutput";

const ORIGIN = "http://treetime.test";

const client = createApiClient({ baseUrl: ORIGIN, fetch: () => Promise.reject(new Error("no request expected")) });

describe("save run output with a desktop host", () => {
  test("the host saves the file and its reply comes back", async () => {
    const saved: unknown[] = [];
    const downloads: string[] = [];

    const host: NonNullable<SaveTarget["host"]> = {
      saveRun: (output) => {
        saved.push(output);

        return Promise.resolve({ kind: "canceled" });
      },
    };

    const reply = await saveRunOutput(
      {
        host,
        client,
        download: (url) => {
          downloads.push(url);
        },
      },
      { id: "r1", path: "a/tree.nwk", name: "tree.nwk" },
    );

    expect([reply, saved, downloads]).toStrictEqual([
      { kind: "canceled" },
      [{ id: "r1", path: "a/tree.nwk", name: "tree.nwk" }],
      [],
    ]);
  });
});

describe("save run output in a browser", () => {
  test("a file starts a download of its file URL", async () => {
    const downloads: Array<[string, string]> = [];

    const reply = await saveRunOutput(
      {
        host: null,
        client,
        download: (url, name) => {
          downloads.push([url, name]);
        },
      },
      { id: "r1", path: "a/tree.nwk", name: "tree.nwk" },
    );

    expect([reply, downloads]).toStrictEqual([
      { kind: "saved" },
      [[`${ORIGIN}/api/runs/r1/file?path=a%2Ftree.nwk`, "tree.nwk"]],
    ]);
  });

  test("an output without a path starts a download of the run archive", async () => {
    const downloads: Array<[string, string]> = [];

    await saveRunOutput(
      {
        host: null,
        client,
        download: (url, name) => {
          downloads.push([url, name]);
        },
      },
      { id: "r1", name: "run.zip" },
    );

    expect(downloads).toStrictEqual([[`${ORIGIN}/api/runs/r1/archive`, "run.zip"]]);
  });
});
