import type { PortReply } from "@neherlab/app-napi";
import { BACKEND_PORT_CHANNEL, BACKEND_PORT_REQUEST_CHANNEL, type BackendStopped } from "@neherlab/app-ui/host";
import { describe, expect, test } from "vitest";

import {
  createPreloadHost,
  invoke,
  windowFetchConnection,
  type IpcRendererLike,
  type WindowLike,
} from "../ipc-renderer";
import type { FetchPort } from "../port-fetch";

describe("ipc_renderer invoke", () => {
  test("a request goes to the prefixed channel and the reply comes back typed", async () => {
    const ipc = fakeIpc({ "treetime:pick-files": ["/data/tree.nwk"] });

    const reply = await invoke(ipc, "pick-files", { title: "Tree", extensions: ["nwk"], multiple: false });

    expect([reply, ipc.invoked]).toStrictEqual([
      ["/data/tree.nwk"],
      [["treetime:pick-files", { title: "Tree", extensions: ["nwk"], multiple: false }]],
    ]);
  });

  test("a reply that fails its schema rejects", async () => {
    const ipc = fakeIpc({ "treetime:pick-folder": 42 });

    await expect(invoke(ipc, "pick-folder", { title: "Runs" })).rejects.toThrow("expected string");
  });

  test("a save reply carries the error response of the back end", async () => {
    const error = { code: "not_found", message: "run `r9` does not exist", causes: [] };
    const ipc = fakeIpc({ "treetime:save-run": { kind: "error", error } });

    await expect(invoke(ipc, "save-run", { id: "r9", name: "r9.zip" })).resolves.toStrictEqual({
      kind: "error",
      error,
    });
  });
});

describe("ipc_renderer preload host", () => {
  test("the host asks for a port on the port request channel", () => {
    const ipc = fakeIpc({});

    createPreloadHost(ipc, { getPathForFile: () => "" }).connectBackend();

    expect(ipc.sent).toStrictEqual([BACKEND_PORT_REQUEST_CHANNEL]);
  });

  test("a stop event reaches the listener after its schema check", () => {
    const ipc = fakeIpc({});
    const stops: BackendStopped[] = [];

    createPreloadHost(ipc, { getPathForFile: () => "" }).onBackendStopped((stop) => {
      stops.push(stop);
    });
    ipc.emit("treetime:backend-stopped", { reason: "crashed", restarts: true });

    expect(stops).toStrictEqual([{ reason: "crashed", restarts: true }]);
  });

  test("the path of a dropped file comes from the file utilities of Electron", () => {
    const host = createPreloadHost(fakeIpc({}), { getPathForFile: (file) => `/home/user/${file.name}` });

    expect(host.pathForFile(new File([], "tree.nwk"))).toBe("/home/user/tree.nwk");
  });
});

describe("ipc_renderer window connection", () => {
  test("the port the preload posts to the window serves the fetch transport, also to listeners added later", () => {
    const target = fakeWindow();
    const host = fakeHost();
    const connection = windowFetchConnection(target, host);
    const port = fakeFetchPort();
    const received: FetchPort[] = [];

    connection.onPort((next) => {
      received.push(next);
    });
    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [port] });
    connection.onPort((next) => {
      received.push(next);
    });

    expect([host.connections, received]).toStrictEqual([1, [port, port]]);
  });

  test("a stopped back end drops its port until the preload posts a new one", () => {
    const target = fakeWindow();
    const host = fakeHost();
    const connection = windowFetchConnection(target, host);
    const received: FetchPort[] = [];

    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [fakeFetchPort()] });
    host.stop({ reason: "crashed", restarts: true });
    connection.onPort((next) => {
      received.push(next);
    });

    expect([host.connections, received]).toStrictEqual([1, []]);
  });

  test("messages from another source, without a port or of another channel are ignored", () => {
    const target = fakeWindow();
    const connection = windowFetchConnection(target, fakeHost());
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
    const host = fakeHost();
    const connection = windowFetchConnection(fakeWindow(), host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    host.stop({ reason: "crashed", restarts: false });

    expect(stops).toStrictEqual([{ reason: "crashed", restarts: false }]);
  });
});

interface FakeHost {
  connections: number;
  connectBackend(): void;
  onBackendStopped(listener: (stop: BackendStopped) => void): void;
  stop(stop: BackendStopped): void;
}

function fakeHost(): FakeHost {
  const stopListeners: Array<(stop: BackendStopped) => void> = [];

  const host: FakeHost = {
    connections: 0,
    connectBackend() {
      host.connections += 1;
    },
    onBackendStopped(listener) {
      stopListeners.push(listener);
    },
    stop(stop) {
      stopListeners.forEach((listener) => {
        listener(stop);
      });
    },
  };

  return host;
}

interface FakeIpc extends IpcRendererLike<PortReply> {
  invoked: unknown[][];
  sent: string[];
  emit(channel: string, payload: unknown): void;
}

function fakeIpc(replies: Record<string, unknown>): FakeIpc {
  const listeners = new Map<string, Array<(event: { ports: PortReply[] }, ...args: unknown[]) => void>>();

  const ipc: FakeIpc = {
    invoked: [],
    sent: [],
    invoke(channel, ...args) {
      ipc.invoked.push([channel, ...args]);

      return Promise.resolve(replies[channel]);
    },
    send(channel) {
      ipc.sent.push(channel);
    },
    on(channel, listener) {
      listeners.set(channel, [...(listeners.get(channel) ?? []), listener]);
    },
    emit(channel, payload) {
      listeners.get(channel)?.forEach((listener) => {
        listener({ ports: [] }, payload);
      });
    },
  };

  return ipc;
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
