import { describe, expect, test } from "vitest";

import { BACKEND_PORT_CHANNEL } from "../channels";
import {
  SaveError,
  createDesktopSaveActions,
  createLocalFiles,
  createWorkspaceShell,
  windowFetchConnection,
  type DesktopShell,
  type WindowLike,
} from "../desktop-shell";
import type { FetchPort } from "../port-fetch";
import type { SaveReply } from "../shell-protocol";

describe("desktop_shell files", () => {
  test("saveRunFile asks the shell to save the file and reports whether it was saved", async () => {
    const shell = fakeShell({ saved: { saved: true } });

    await expect(createDesktopSaveActions(shell).saveRunFile("r1", "out/a.nwk", "a.nwk")).resolves.toBe(true);
    expect(shell.saves).toStrictEqual([{ id: "r1", path: "out/a.nwk", name: "a.nwk" }]);
  });

  test("saveRunArchive resolves false when the user cancels the save dialog", async () => {
    const shell = fakeShell({ saved: { saved: false } });

    await expect(createDesktopSaveActions(shell).saveRunArchive("r1", "run.zip")).resolves.toBe(false);
    expect(shell.saves).toStrictEqual([{ id: "r1", name: "run.zip" }]);
  });

  test("a save the back end refuses rejects with its typed error", async () => {
    const response = { code: "invalid_request", message: "file path `../x` must name a file", causes: [] };
    const shell = fakeShell({ saved: { error: JSON.stringify(response) } });

    const error: unknown = await createDesktopSaveActions(shell)
      .saveRunFile("r1", "../x", "x")
      .catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(SaveError);
    expect(error).toMatchObject({ response, message: "file path `../x` must name a file" });
  });

  test("an untyped save error rejects as an internal error with its message", async () => {
    const shell = fakeShell({ saved: { error: "permission denied" } });

    await expect(createDesktopSaveActions(shell).saveRunArchive("r1", "run.zip")).rejects.toMatchObject({
      response: { code: "internal_error", message: "permission denied", causes: [] },
    });
  });

  test("pickFiles validates the paths the shell returns", async () => {
    const files = createLocalFiles(fakeShell({ picked: ["/data/tree.nwk"] }));

    await expect(files.pickFiles({ title: "Tree", extensions: ["nwk"], multiple: false })).resolves.toStrictEqual([
      "/data/tree.nwk",
    ]);
  });

  test("a malformed pick result rejects", async () => {
    const files = createLocalFiles(fakeShell({ picked: [1] }));

    await expect(files.pickFiles({ title: "Tree", extensions: [], multiple: false })).rejects.toThrow(
      "expected string",
    );
  });
});

describe("desktop_shell workspace", () => {
  test("pickFolder returns the folder the user chose", async () => {
    const workspace = createWorkspaceShell(fakeShell({ folder: "/data/runs" }));

    await expect(workspace.pickFolder({ title: "Runs folder" })).resolves.toBe("/data/runs");
  });

  test("pickFolder returns null when the user cancels the dialog", async () => {
    const workspace = createWorkspaceShell(fakeShell({ folder: null }));

    await expect(workspace.pickFolder({ title: "Runs folder" })).resolves.toBeNull();
  });

  test("a malformed folder result rejects", async () => {
    const workspace = createWorkspaceShell(fakeShell({ folder: ["/data/runs"] }));

    await expect(workspace.pickFolder({ title: "Runs folder" })).rejects.toThrow("expected string");
  });

  test("restartBackend asks the shell to restart the back end once", async () => {
    const shell = fakeShell({});

    await createWorkspaceShell(shell).restartBackend();

    expect(shell.restarts).toBe(1);
  });
});

describe("desktop_shell window connection", () => {
  test("the port the preload posts to the window serves the fetch transport, also to listeners added later", () => {
    const target = fakeWindow();
    const shell = fakeShell({});
    const connection = windowFetchConnection(target, shell);
    const port = fakeFetchPort();
    const received: FetchPort[] = [];

    connection.onPort((next) => {
      received.push(next);
    });
    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [port] });
    connection.onPort((next) => {
      received.push(next);
    });

    expect([shell.connections, received]).toStrictEqual([1, [port, port]]);
  });

  test("a stopped back end drops its port until the preload posts a new one", () => {
    const target = fakeWindow();
    const shell = fakeShell({});
    const connection = windowFetchConnection(target, shell);
    const received: FetchPort[] = [];

    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [fakeFetchPort()] });
    shell.stop("crashed", true);
    connection.onPort((next) => {
      received.push(next);
    });

    expect([shell.connections, received]).toStrictEqual([1, []]);
  });

  test("messages from another source, without a port or of another channel are ignored", () => {
    const target = fakeWindow();
    const connection = windowFetchConnection(target, fakeShell({}));
    const received: FetchPort[] = [];

    connection.onPort((next) => {
      received.push(next);
    });
    target.emit({ source: {}, data: { channel: BACKEND_PORT_CHANNEL }, ports: [fakeFetchPort()] });
    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [] });
    target.emit({ source: target, data: "treetime", ports: [fakeFetchPort()] });
    target.emit({ source: target, data: { channel: "treetime:other" }, ports: [fakeFetchPort()] });

    expect(received).toStrictEqual([]);
  });

  test("a stop of the back end reaches the listeners with its reason", () => {
    const shell = fakeShell({});
    const connection = windowFetchConnection(fakeWindow(), shell);
    const stops: Array<[string, boolean]> = [];

    connection.onStopped((reason, restarts) => {
      stops.push([reason, restarts]);
    });
    shell.stop("crashed", false);

    expect(stops).toStrictEqual([["crashed", false]]);
  });
});

interface FakeShell extends DesktopShell {
  connections: number;
  restarts: number;
  saves: unknown[];
  stop(reason: string, restarts: boolean): void;
}

function fakeShell({
  picked = [],
  folder = null,
  saved = { saved: false },
}: {
  picked?: unknown;
  folder?: unknown;
  saved?: SaveReply;
}): FakeShell {
  const stopListeners: Array<(reason: string, restarts: boolean) => void> = [];

  const shell: FakeShell = {
    connections: 0,
    restarts: 0,
    saves: [],
    connectBackend() {
      shell.connections += 1;
    },
    onBackendStopped(listener) {
      stopListeners.push(listener);
    },
    stop(reason, restarts) {
      stopListeners.forEach((listener) => {
        listener(reason, restarts);
      });
    },
    pickFiles: () => Promise.resolve(picked),
    pickFolder: () => Promise.resolve(folder),
    restartBackend: () => {
      shell.restarts += 1;

      return Promise.resolve();
    },
    pathForFile: () => "",
    saveRunFile: (request) => {
      shell.saves.push(request);

      return Promise.resolve(saved);
    },
    saveRunArchive: (request) => {
      shell.saves.push(request);

      return Promise.resolve(saved);
    },
  };

  return shell;
}

type WindowMessage = Parameters<Parameters<WindowLike["addEventListener"]>[1]>[0];

function fakeWindow() {
  const listeners: Array<(event: WindowMessage) => void> = [];

  return {
    addEventListener(_type: "message", listener: (event: WindowMessage) => void) {
      listeners.push(listener);
    },
    emit(event: WindowMessage) {
      listeners.forEach((listener) => {
        listener(event);
      });
    },
  };
}

function fakeFetchPort(): FetchPort {
  return { postMessage: () => undefined, addEventListener: () => undefined, start: () => undefined };
}
