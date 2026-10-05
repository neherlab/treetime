import type { ErrorResponse } from "@neherlab/app-contracts";
import { ApiError } from "@neherlab/app-contracts/client";
import type { PortMessage, PortReply } from "@neherlab/app-napi";
import type { BackendStopped } from "@neherlab/app-ui/host";
import { describe, expect, test, vi } from "vitest";

import { createPortFetch, type FetchConnection, type FetchPort } from "../port-fetch";

const ORIGIN = "http://treetime.desktop";

describe("port_fetch requests", () => {
  test("a request made before the port connects is sent once it connects", async () => {
    const backend = fakeBackend();
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    const before = [...backend.posted];

    backend.connect();
    backend.reply({ kind: "end", seq: 0 });

    await expect(response).rejects.toThrow(new TypeError("the back end ended the exchange without a response head"));
    expect([before, backend.posted]).toStrictEqual([
      [],
      [{ kind: "request", request: { seq: 0, method: "GET", url: "/api/version", headers: [] } }],
    ]);
  });

  test("a request sends its method, path, query, headers and body", async () => {
    const backend = fakeBackend({ connected: true });

    void createPortFetch(backend.connection)(`${ORIGIN}/api/runs?from=2`, {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: '{"command":"clock"}',
    });
    await sent(backend, 1);

    expect(backend.posted).toStrictEqual([
      {
        kind: "request",
        request: {
          seq: 0,
          method: "POST",
          url: "/api/runs?from=2",
          headers: [{ name: "content-type", value: "application/json" }],
          body: new TextEncoder().encode('{"command":"clock"}'),
        },
      },
    ]);
  });

  test("a binary body arrives byte for byte, also when it is not valid UTF-8", async () => {
    const backend = fakeBackend({ connected: true });
    const bytes = new Uint8Array([0xff, 0xfe, 0x00, 0xc3, 0x28, 0x80]);

    void createPortFetch(backend.connection)(`${ORIGIN}/api/runs/r1/inputs/blob.bin`, { method: "PUT", body: bytes });
    await sent(backend, 1);

    expect(backend.posted.map((message) => message.kind === "request" && message.request.body)).toStrictEqual([bytes]);
  });

  test("a request without a body sends no body", async () => {
    const backend = fakeBackend({ connected: true });

    void createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    await vi.waitFor(() => {
      expect(backend.posted).toHaveLength(1);
    });

    expect(backend.posted).toStrictEqual([
      { kind: "request", request: { seq: 0, method: "GET", url: "/api/version", headers: [] } },
    ]);
  });
});

describe("port_fetch responses", () => {
  test("the head resolves the response and the chunks stream into its body without decoding", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    await sent(backend, 1);

    backend.reply({ kind: "head", seq: 0, status: 201, headers: [{ name: "content-type", value: "text/plain" }] });
    const resolved = await response;
    const reader = resolved.body?.getReader();
    backend.reply({ kind: "chunk", seq: 0, data: new Uint8Array([0x61, 0xe2, 0x82]) });
    const first = await reader?.read();
    backend.reply({ kind: "chunk", seq: 0, data: new Uint8Array([0xac]) });
    backend.reply({ kind: "end", seq: 0 });
    const second = await reader?.read();
    const done = await reader?.read();

    expect([resolved.status, resolved.headers.get("content-type"), first, second, done?.done]).toStrictEqual([
      201,
      "text/plain",
      { done: false, value: new Uint8Array([0x61, 0xe2, 0x82]) },
      { done: false, value: new Uint8Array([0xac]) },
      true,
    ]);
  });

  test("a complete body reads as the text of its chunks", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    await sent(backend, 1);

    backend.reply({ kind: "head", seq: 0, status: 200, headers: [] });
    backend.reply({ kind: "chunk", seq: 0, data: new TextEncoder().encode('{"version":') });
    backend.reply({ kind: "chunk", seq: 0, data: new TextEncoder().encode('"2.0.0"}') });
    backend.reply({ kind: "end", seq: 0 });

    await expect((await response).json()).resolves.toStrictEqual({ version: "2.0.0" });
  });

  test("a 204 answer resolves with an empty body", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/runs/r1`, { method: "DELETE" });
    await sent(backend, 1);

    backend.reply({ kind: "head", seq: 0, status: 204, headers: [] });
    backend.reply({ kind: "end", seq: 0 });
    const resolved = await response;

    expect([resolved.status, resolved.body]).toStrictEqual([204, null]);
  });

  test("a reset before the head rejects with its message", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    await sent(backend, 1);

    backend.reply({ kind: "reset", seq: 0, message: "the back end stopped the operation" });

    await expect(response).rejects.toThrow(new TypeError("the back end stopped the operation"));
  });

  test("a reset after the head errors the body, so a truncated body never looks complete", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`);
    await sent(backend, 1);

    backend.reply({ kind: "head", seq: 0, status: 200, headers: [] });
    const resolved = await response;
    backend.reply({ kind: "chunk", seq: 0, data: new Uint8Array([104]) });
    backend.reply({ kind: "reset", seq: 0, message: "When reading the response body: disk gone" });

    await expect(resolved.text()).rejects.toThrow("When reading the response body: disk gone");
  });
});

describe("port_fetch aborts", () => {
  test("an abort before the head rejects with the abort reason and aborts the exchange", async () => {
    const backend = fakeBackend({ connected: true });
    const controller = new AbortController();
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`, { signal: controller.signal });
    await sent(backend, 1);

    controller.abort(new Error("the view closed"));

    await expect(response).rejects.toThrow("the view closed");
    expect(backend.posted.at(-1)).toStrictEqual({ kind: "abort", seq: 0 });
  });

  test("an abort while the body streams errors the body and aborts the exchange", async () => {
    const backend = fakeBackend({ connected: true });
    const controller = new AbortController();

    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/runs/r1/events`, {
      signal: controller.signal,
    });

    await sent(backend, 1);
    backend.reply({ kind: "head", seq: 0, status: 200, headers: [] });
    const resolved = await response;

    controller.abort(new Error("the view closed"));

    await expect(resolved.text()).rejects.toThrow("the view closed");
    expect(backend.posted.at(-1)).toStrictEqual({ kind: "abort", seq: 0 });
  });

  test("cancelling the body aborts the exchange", async () => {
    const backend = fakeBackend({ connected: true });
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/runs/r1/events`);
    await sent(backend, 1);
    backend.reply({ kind: "head", seq: 0, status: 200, headers: [] });

    await (await response).body?.cancel();

    expect(backend.posted.at(-1)).toStrictEqual({ kind: "abort", seq: 0 });
  });

  test("an abort of a request that waits for the port sends nothing", async () => {
    const backend = fakeBackend();
    const controller = new AbortController();
    const response = createPortFetch(backend.connection)(`${ORIGIN}/api/version`, { signal: controller.signal });

    controller.abort(new Error("the view closed"));
    backend.connect();

    await expect(response).rejects.toThrow("the view closed");
    expect(backend.posted).toStrictEqual([]);
  });

  test("an aborted signal rejects before anything is sent", async () => {
    const backend = fakeBackend({ connected: true });

    await expect(
      createPortFetch(backend.connection)(`${ORIGIN}/api/version`, { signal: AbortSignal.abort(new Error("late")) }),
    ).rejects.toThrow("late");
    expect(backend.posted).toStrictEqual([]);
  });
});

describe("port_fetch back end stops", () => {
  test("a restart rejects pending requests and errors open bodies; later requests wait for the new port", async () => {
    const backend = fakeBackend({ connected: true });
    const portFetch = createPortFetch(backend.connection);
    const pending = portFetch(`${ORIGIN}/api/version`);
    const streaming = portFetch(`${ORIGIN}/api/runs/r1/events`);
    await sent(backend, 2);
    backend.reply({ kind: "head", seq: 1, status: 200, headers: [] });
    const open = await streaming;

    backend.stop({ reason: "the back end stopped with exit code 1 and restarts", restarts: true });
    await expect(pending).rejects.toThrow(new TypeError("the back end stopped with exit code 1 and restarts"));
    await expect(open.text()).rejects.toThrow("the back end stopped with exit code 1 and restarts");
    const later = portFetch(`${ORIGIN}/api/version`);
    const postedWhileStopped = backend.posted.length;
    backend.connect();
    await sent(backend, 3);
    backend.reply({ kind: "head", seq: 2, status: 200, headers: [] });
    backend.reply({ kind: "end", seq: 2 });

    await expect(later.then((response) => response.status)).resolves.toBe(200);
    expect([postedWhileStopped, backend.posted.at(-1)]).toStrictEqual([
      2,
      { kind: "request", request: { seq: 2, method: "GET", url: "/api/version", headers: [] } },
    ]);
  });

  test("a back end that does not restart rejects every later request with its reason", async () => {
    const backend = fakeBackend({ connected: true });
    const portFetch = createPortFetch(backend.connection);

    backend.stop({ reason: "the back end stopped too often", restarts: false });

    await expect(portFetch(`${ORIGIN}/api/version`)).rejects.toThrow("the back end stopped too often");
    expect(backend.posted).toStrictEqual([]);
  });

  test("a back end that could not start rejects pending and later requests with its error response", async () => {
    const backend = fakeBackend();
    const portFetch = createPortFetch(backend.connection);
    const pending = portFetch(`${ORIGIN}/api/version`);

    const error: ErrorResponse = {
      code: "internal_error",
      message: "When opening the runs folder '/r'",
      causes: ["denied"],
    };

    backend.stop({ reason: "the back end could not start", restarts: false, error });

    const expected = new ApiError(500, error);

    const rejections = await Promise.all(
      [pending, portFetch(`${ORIGIN}/api/runs`)].map(async (request) => request.catch((reason: unknown) => reason)),
    );

    expect(rejections).toStrictEqual([expected, expected]);
    expect(rejections.map((rejection) => rejection instanceof ApiError && rejection.status)).toStrictEqual([500, 500]);
  });

  test("a new port after a stop without restart accepts requests again", async () => {
    const backend = fakeBackend({ connected: true });
    const portFetch = createPortFetch(backend.connection);

    backend.stop({ reason: "the back end stopped too often", restarts: false });
    backend.connect();
    void portFetch(`${ORIGIN}/api/version`);
    await sent(backend, 1);

    expect(backend.posted).toStrictEqual([
      { kind: "request", request: { seq: 0, method: "GET", url: "/api/version", headers: [] } },
    ]);
  });
});

interface FakeBackend {
  connection: FetchConnection;
  posted: PortMessage[];
  connect(): void;
  reply(reply: PortReply): void;
  stop(stop: BackendStopped): void;
}

function fakeBackend({ connected = false } = {}): FakeBackend {
  const portListeners: Array<(port: FetchPort) => void> = [];
  const stoppedListeners: Array<(stop: BackendStopped) => void> = [];
  let deliver: (reply: PortReply) => void = () => undefined;

  const backend: FakeBackend = {
    posted: [],
    connection: {
      onPort(listener) {
        portListeners.push(listener);

        if (connected) {
          backend.connect();
        }
      },
      onStopped(listener) {
        stoppedListeners.push(listener);
      },
    },
    connect() {
      const port: FetchPort = {
        postMessage(message) {
          backend.posted.push(message);
        },
        addEventListener(_type, listener) {
          deliver = (data) => {
            listener({ data });
          };
        },
        start: () => undefined,
      };

      portListeners.forEach((listener) => {
        listener(port);
      });
    },
    reply(reply) {
      deliver(reply);
    },
    stop(stop) {
      deliver = () => undefined;
      stoppedListeners.forEach((listener) => {
        listener(stop);
      });
    },
  };

  return backend;
}

async function sent(backend: FakeBackend, count: number): Promise<void> {
  await vi.waitFor(() => {
    expect(backend.posted).toHaveLength(count);
  });
}
