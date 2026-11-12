import { createApiClient } from "@neherlab/app-contracts/client";

import { browserFilePicker, type ObjectUrls, saveRunOutput, type SaveTarget } from "../saveRunOutput";

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
        pickFile: null,
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
  test("a file downloads its content through an object URL under the file name", async () => {
    const server = fakeServer({ "/api/runs/r1/file?path=a%2Ftree.nwk": "(A,B);" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const downloads: Download[] = [];

    const reply = await saveRunOutput(
      {
        host: null,
        client: server.client,
        pickFile: null,
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

  test("an output without a path downloads the run package", async () => {
    const server = fakeServer({ "/api/runs/r1/package": "PK" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const downloads: Download[] = [];

    await saveRunOutput(
      {
        host: null,
        client: server.client,
        pickFile: null,
        download: (url, name) => {
          downloads.push({ url, name });
        },
        objectUrls,
        setTimer: timers.set,
      },
      { id: "r1", name: "run.zip" },
    );

    expect([server.requested, await objectUrls.contents(downloads)]).toStrictEqual([
      ["/api/runs/r1/package"],
      [{ url: "blob:treetime/0", name: "run.zip", content: "PK" }],
    ]);
  });

  test("the object URL stays valid after the click and is revoked when its timer fires", async () => {
    const server = fakeServer({ "/api/runs/r1/package": "PK" });
    const objectUrls = new FakeObjectUrls();
    const timers = new FakeTimers();
    const urls: string[] = [];

    await saveRunOutput(
      {
        host: null,
        client: server.client,
        pickFile: null,
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
        client: server.client,
        pickFile: null,
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

describe("save run output in a browser with a save file picker", () => {
  test("the output streams into the picked file without an object URL", async () => {
    const server = fakeServer({ "/api/runs/r1/package": "PK" });
    const objectUrls = new FakeObjectUrls();
    const files = new FakeFiles();

    const reply = await saveRunOutput(
      { host: null, client: server.client, pickFile: files.pick, download: noDownload, objectUrls, setTimer: noTimer },
      { id: "r1", name: "ebola-1a2b3c4d.zip" },
    );

    expect([reply, server.requested, files.saved(), objectUrls.liveCount()]).toStrictEqual([
      { kind: "saved" },
      ["/api/runs/r1/package"],
      [{ name: "ebola-1a2b3c4d.zip", content: "PK", state: "closed" }],
      0,
    ]);
  });

  test("a canceled picker requests nothing", async () => {
    const server = fakeServer({ "/api/runs/r1/package": "PK" });
    const objectUrls = new FakeObjectUrls();

    const reply = await saveRunOutput(
      {
        host: null,
        client: server.client,
        pickFile: () => Promise.resolve(undefined),
        download: noDownload,
        objectUrls,
        setTimer: noTimer,
      },
      { id: "r1", name: "run.zip" },
    );

    expect([reply, server.requested]).toStrictEqual([{ kind: "canceled" }, []]);
  });

  test("a file the back end cannot send aborts the picked file and fails with its error", async () => {
    const server = fakeServer({});
    const objectUrls = new FakeObjectUrls();
    const files = new FakeFiles();

    const save = saveRunOutput(
      { host: null, client: server.client, pickFile: files.pick, download: noDownload, objectUrls, setTimer: noTimer },
      { id: "r1", path: "missing.nwk", name: "missing.nwk" },
    );

    await expect(save).rejects.toThrow("no such file");
    expect(files.saved()).toStrictEqual([{ name: "missing.nwk", content: "", state: "aborted" }]);
  });
});

describe("browser file picker", () => {
  test("a browser without a save file picker or outside a secure context has none", () => {
    const showSaveFilePicker = () => Promise.reject(new Error("no picker expected"));

    expect([
      browserFilePicker({ isSecureContext: true }),
      browserFilePicker({ isSecureContext: false, showSaveFilePicker }),
    ]).toStrictEqual([null, null]);
  });

  test("a dismissed picker counts as canceled", async () => {
    const pick = browserFilePicker({
      isSecureContext: true,
      showSaveFilePicker: () => Promise.reject(new DOMException("The user aborted a request.", "AbortError")),
    });

    expect(await pick?.("run.zip")).toBeUndefined();
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

function noDownload(): void {
  throw new Error("no download expected");
}

interface SavedFile {
  name: string;
  content: string;
  state: "open" | "closed" | "aborted";
}

interface PickedFile {
  file: SavedFile;
  chunks: Uint8Array[];
}

class FakeFiles {
  readonly #files: PickedFile[] = [];

  readonly pick = (name: string): Promise<WritableStream<Uint8Array>> => {
    const entry: PickedFile = { file: { name, content: "", state: "open" }, chunks: [] };

    this.#files.push(entry);

    return Promise.resolve(
      new WritableStream<Uint8Array>({
        write: (chunk) => {
          entry.chunks.push(chunk);
        },
        close: () => {
          entry.file.state = "closed";
        },
        abort: () => {
          entry.file.state = "aborted";
        },
      }),
    );
  };

  saved(): SavedFile[] {
    const decoder = new TextDecoder();

    return this.#files.map(({ file, chunks }) => ({
      ...file,
      content: chunks.map((chunk) => decoder.decode(chunk, { stream: true })).join("") + decoder.decode(),
    }));
  }
}
