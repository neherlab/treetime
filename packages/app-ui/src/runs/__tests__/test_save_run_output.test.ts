import { createApiClient } from "@neherlab/app-contracts/client";

import { TAB_BACKEND } from "../../browser-backend";
import { type ObjectUrls, saveRunOutput, type SaveTarget } from "../saveRunOutput";

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
        browserBackend: null,
        client,
        objectUrls: URL,
        setTimer: noTimer,
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
        browserBackend: null,
        client,
        objectUrls: URL,
        setTimer: noTimer,
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
        browserBackend: null,
        client,
        objectUrls: URL,
        setTimer: noTimer,
        download: (url, name) => {
          downloads.push([url, name]);
        },
      },
      { id: "r1", name: "run.zip" },
    );

    expect(downloads).toStrictEqual([[`${ORIGIN}/api/runs/r1/archive`, "run.zip"]]);
  });
});

describe("save run output with a browser back end", () => {
  test("a file downloads its content through an object URL under the file name", async () => {
    const server = fakeServer({ "/api/runs/r1/file?path=a%2Ftree.nwk": "(A,B);" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const downloads: Download[] = [];

    const reply = await saveRunOutput(
      {
        host: null,
        browserBackend: TAB_BACKEND,
        client: server.client,
        download: (url, name) => {
          downloads.push({ url, name });
        },
        objectUrls,
        setTimer: timers.set,
      },
      { id: "r1", path: "a/tree.nwk", name: "tree.nwk" },
    );

    expect([reply, server.requested, await objectUrls.contents(downloads)]).toStrictEqual([
      { kind: "saved" },
      ["/api/runs/r1/file?path=a%2Ftree.nwk"],
      [{ url: "blob:treetime/0", name: "tree.nwk", content: "(A,B);" }],
    ]);
  });

  test("an output without a path downloads the run archive", async () => {
    const server = fakeServer({ "/api/runs/r1/archive": "PK" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const downloads: Download[] = [];

    await saveRunOutput(
      {
        host: null,
        browserBackend: TAB_BACKEND,
        client: server.client,
        download: (url, name) => {
          downloads.push({ url, name });
        },
        objectUrls,
        setTimer: timers.set,
      },
      { id: "r1", name: "run.zip" },
    );

    expect([server.requested, await objectUrls.contents(downloads)]).toStrictEqual([
      ["/api/runs/r1/archive"],
      [{ url: "blob:treetime/0", name: "run.zip", content: "PK" }],
    ]);
  });

  test("the object URL stays valid after the click and is revoked when its timer fires", async () => {
    const server = fakeServer({ "/api/runs/r1/archive": "PK" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const urls: string[] = [];

    await saveRunOutput(
      {
        host: null,
        browserBackend: TAB_BACKEND,
        client: server.client,
        download: (url) => {
          urls.push(url);
        },
        objectUrls,
        setTimer: timers.set,
      },
      { id: "r1", name: "run.zip" },
    );
    const liveAfterClick = urls.map((url) => objectUrls.isLive(url));
    timers.fireAll();

    expect([liveAfterClick, urls.map((url) => objectUrls.isLive(url))]).toStrictEqual([[true], [false]]);
  });

  test("a file the back end cannot send starts no download and fails with its error", async () => {
    const server = fakeServer({});
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const downloads: string[] = [];

    const save = saveRunOutput(
      {
        host: null,
        browserBackend: TAB_BACKEND,
        client: server.client,
        download: (url) => {
          downloads.push(url);
        },
        objectUrls,
        setTimer: timers.set,
      },
      { id: "r1", path: "missing.nwk", name: "missing.nwk" },
    );

    await expect(save).rejects.toThrow("no such file");
    expect([downloads, objectUrls.liveCount()]).toStrictEqual([[], 0]);
  });
});

interface FakeServer {
  client: SaveTarget["client"];
  requested: string[];
}

function fakeServer(files: Record<string, string>): FakeServer {
  const requested: string[] = [];

  return {
    requested,
    client: createApiClient({
      baseUrl: ORIGIN,
      fetch: (input) => {
        const url = new URL(input instanceof Request ? input.url : String(input));
        const path = `${url.pathname}${url.search}`;
        requested.push(path);
        const content = files[path];

        return Promise.resolve(
          content === undefined
            ? Response.json({ message: "no such file", causes: [] }, { status: 404 })
            : new Response(content, { headers: { "content-type": "application/octet-stream" } }),
        );
      },
    }),
  };
}

interface Download {
  url: string;
  name: string;
}

class FakeObjectUrls implements ObjectUrls {
  readonly #created = new Map<string, Blob>();
  readonly #live = new Set<string>();

  createObjectURL(blob: Blob): string {
    const url = `blob:treetime/${this.#created.size}`;

    this.#created.set(url, blob);
    this.#live.add(url);

    return url;
  }

  revokeObjectURL(url: string): void {
    this.#live.delete(url);
  }

  isLive(url: string): boolean {
    return this.#live.has(url);
  }

  liveCount(): number {
    return this.#live.size;
  }

  contents(downloads: Download[]): Promise<Array<Download & { content: string | undefined }>> {
    return Promise.all(
      downloads.map(async ({ url, name }) => ({ url, name, content: await this.#created.get(url)?.text() })),
    );
  }
}

class FakeTimers {
  readonly #pending: Array<() => void> = [];

  readonly set = (callback: () => void): void => {
    this.#pending.push(callback);
  };

  fireAll(): void {
    for (const callback of this.#pending.splice(0)) {
      callback();
    }
  }
}

function noTimer(): void {
  throw new Error("no timer expected");
}
