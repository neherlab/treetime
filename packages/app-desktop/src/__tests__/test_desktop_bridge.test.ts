import { BridgeError, CancelledError, RunEndedError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import type { BackendReply, BackendRequest, ClientEndpoint, PortLike } from "../backend-protocol";
import { BACKEND_PORT_CHANNEL, BACKEND_STOPPED_CHANNEL } from "../channels";
import {
  createDesktopBridge,
  createLocalFiles,
  windowBackendConnection,
  type BackendConnection,
  type DesktopShell,
  type WindowMessage,
} from "../desktop-bridge";

const RECORD = {
  id: "r1",
  title: "clock",
  command: "clock",
  config: { tree: "t.nwk" },
  status: "running",
  pinned: false,
  created_at: "2026-09-25T10:00:00Z",
  started_at: null,
  finished_at: null,
  duration_seconds: null,
  treetime_version: "1.0.0",
  inputs: [],
  config_hash: null,
  changed_settings: [],
  headline: {},
  output_files: [],
  error: null,
};

const OUTCOME = { command: "clock", output_files: [{ path: "/runs/r1/out/clock.nwk", kind: "nwk" }] };

describe("desktop_bridge operations", () => {
  test("version sends the version operation and validates the result", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "result", seq: request.seq, json: JSON.stringify({ version: "2.0.0" }) });
    });

    await expect(createDesktopBridge(backend.connection, fakeShell({ picked: [] })).version()).resolves.toStrictEqual({
      version: "2.0.0",
    });
    expect(backend.requests).toStrictEqual([
      { kind: "call", seq: 0, request: JSON.stringify({ operation: "version", args: {} }) },
    ]);
  });

  test("checkConfig sends the check-config operation with its request", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({
        kind: "result",
        seq: request.seq,
        json: JSON.stringify({
          status: "valid",
          command: "clock",
          config: {},
          code: { command_line: [], command_line_text: "", yaml: [], yaml_text: "" },
          checks: [],
        }),
      });
    });

    await createDesktopBridge(backend.connection, fakeShell({ picked: [] })).checkConfig({
      command: "clock",
      text: "tree: t",
    });
    expect(backend.calls()).toStrictEqual([
      { operation: "check-config", args: { request: { command: "clock", text: "tree: t" } } },
    ]);
  });

  test("startRun sends the start-run operation with the replacement configuration", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "result", seq: request.seq, json: JSON.stringify(RECORD) });
    });

    await createDesktopBridge(backend.connection, fakeShell({ picked: [] })).startRun("r1", {
      config: { tree: "/data/t.nwk" },
    });
    expect(backend.calls()).toStrictEqual([
      { operation: "start-run", args: { id: "r1", request: { config: { tree: "/data/t.nwk" } } } },
    ]);
  });

  test("a malformed result rejects with the validation error of the bridge", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "result", seq: request.seq, json: JSON.stringify({ version: 2 }) });
    });

    await expect(createDesktopBridge(backend.connection, fakeShell({ picked: [] })).version()).rejects.toThrow(
      "expected string",
    );
  });

  test("calls made before the back end connects are sent once it connects", async () => {
    const backend = fakeBackend(
      (request, reply) => {
        reply({ kind: "result", seq: request.seq, json: JSON.stringify({ version: "2.0.0" }) });
      },
      { connected: false },
    );

    const version = createDesktopBridge(backend.connection, fakeShell({ picked: [] })).version();
    await Promise.resolve();
    expect(backend.requests).toStrictEqual([]);

    backend.connect();

    await expect(version).resolves.toStrictEqual({ version: "2.0.0" });
  });

  test("uploadInput rejects because the desktop app reads local paths", async () => {
    const backend = fakeBackend(() => undefined);
    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).uploadInput("r1", "t.nwk", new Blob(["x"])),
    ).rejects.toThrow(
      "the desktop application reads inputs from local file paths; name the files in the run configuration",
    );
  });
});

describe("desktop_bridge errors", () => {
  test("a typed error of the back end rejects with a BridgeError carrying its class and causes", async () => {
    const response = { code: "not_found", message: "When reading run `r9`", causes: ["no run with id `r9`"] };

    const backend = fakeBackend((request, reply) => {
      reply({ kind: "error", seq: request.seq, error: JSON.stringify(response) });
    });

    const error = await createDesktopBridge(backend.connection, fakeShell({ picked: [] }))
      .getRun("r9")
      .catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(BridgeError);
    expect(error).toMatchObject({ message: "When reading run `r9`: no run with id `r9`", response });
  });

  test("an untyped error rejects as an internal error with its message", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "error", seq: request.seq, error: "addon failed to load" });
    });

    await expect(createDesktopBridge(backend.connection, fakeShell({ picked: [] })).version()).rejects.toMatchObject({
      response: { code: "internal_error", message: "addon failed to load", causes: [] },
    });
  });

  test("a stopped back end rejects the calls it did not answer", async () => {
    const backend = fakeBackend(() => undefined);
    const version = createDesktopBridge(backend.connection, fakeShell({ picked: [] })).version();
    await Promise.resolve();

    backend.stop("the back end stopped with exit code 134 and restarts; the request was not answered");

    await expect(version).rejects.toMatchObject({
      response: {
        code: "internal_error",
        message: "the back end stopped with exit code 134 and restarts; the request was not answered",
        causes: [],
      },
    });
  });
});

describe("desktop_bridge run events", () => {
  test("followRun subscribes, forwards the events and unsubscribes at the terminal event", async () => {
    const backend = fakeBackend((request, reply) => {
      if (request.kind === "subscribe") {
        reply({ kind: "event", seq: request.seq, json: runEvent(3, "log", { level: "info", message: "hello" }) });
        reply({
          kind: "event",
          seq: request.seq,
          json: runEvent(4, "terminal", { status: "cancelled", job_id: "r1" }),
        });
      }
    });

    const seen: number[] = [];

    const terminal = await createDesktopBridge(backend.connection, fakeShell({ picked: [] })).followRun("r1", {
      from: 3,
      onEvent: (event) => {
        seen.push(event.seq);
      },
    });

    expect(terminal).toStrictEqual({ status: "cancelled", job_id: "r1" });
    expect(seen).toStrictEqual([3, 4]);
    expect(backend.requests).toStrictEqual([
      { kind: "subscribe", seq: 0, id: "r1", from: 3 },
      { kind: "unsubscribe", seq: 0 },
    ]);
  });

  test("aborting followRun unsubscribes the run events", async () => {
    const backend = fakeBackend(() => undefined);
    const controller = new AbortController();

    const following = createDesktopBridge(backend.connection, fakeShell({ picked: [] })).followRun("r1", {
      signal: controller.signal,
    });

    await Promise.resolve();
    controller.abort();

    await expect(following).rejects.toBeInstanceOf(RunEndedError);
    expect(backend.requests).toStrictEqual([
      { kind: "subscribe", seq: 0, id: "r1", from: 0 },
      { kind: "unsubscribe", seq: 0 },
    ]);
  });

  test("a restarted back end resumes the run events after the last event received", async () => {
    let attempt = 0;

    const backend = fakeBackend((request, reply) => {
      if (request.kind !== "subscribe") {
        return;
      }

      attempt += 1;

      if (attempt === 1) {
        reply({ kind: "event", seq: request.seq, json: runEvent(0, "started", { job_id: "r1", command: "clock" }) });
      } else {
        reply({
          kind: "event",
          seq: request.seq,
          json: runEvent(1, "terminal", { status: "interrupted", job_id: "r1" }),
        });
      }
    });

    const following = createDesktopBridge(backend.connection, fakeShell({ picked: [] })).followRun("r1");
    await Promise.resolve();

    backend.stop("the back end stopped with exit code 134 and restarts; the request was not answered");
    backend.connect();

    await expect(following).resolves.toStrictEqual({ status: "interrupted", job_id: "r1" });
    expect(backend.requests.filter((request) => request.kind === "subscribe")).toStrictEqual([
      { kind: "subscribe", seq: 0, id: "r1", from: 0 },
      { kind: "subscribe", seq: 0, id: "r1", from: 1 },
    ]);
  });

  test("a command creates a run and resolves with the outcome", async () => {
    const backend = fakeBackend((request, reply) => {
      if (request.kind === "call") {
        reply({ kind: "result", seq: request.seq, json: JSON.stringify(RECORD) });
      } else if (request.kind === "subscribe") {
        reply({ kind: "event", seq: request.seq, json: runEvent(0, "started", { job_id: "r1", command: "clock" }) });
        reply({
          kind: "event",
          seq: request.seq,
          json: runEvent(1, "terminal", { status: "ok", job_id: "r1", result: OUTCOME }),
        });
      }
    });

    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).clock({ tree: "t.nwk" }),
    ).resolves.toStrictEqual(OUTCOME);
  });

  test("a cancelled run rejects the command with a CancelledError", async () => {
    const backend = fakeBackend((request, reply) => {
      if (request.kind === "call") {
        reply({ kind: "result", seq: request.seq, json: JSON.stringify(RECORD) });
      } else if (request.kind === "subscribe") {
        reply({
          kind: "event",
          seq: request.seq,
          json: runEvent(0, "terminal", { status: "cancelled", job_id: "r1" }),
        });
      }
    });

    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).clock({ tree: "t.nwk" }),
    ).rejects.toBeInstanceOf(CancelledError);
  });

  test("a failed subscription rejects with the error of the back end", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({
        kind: "error",
        seq: request.seq,
        error: JSON.stringify({ code: "not_found", message: "no run with id `r9`", causes: [] }),
      });
    });

    await expect(createDesktopBridge(backend.connection, fakeShell({ picked: [] })).followRun("r9")).rejects.toThrow(
      "no run with id `r9`",
    );
  });

  test("a malformed run event rejects with the validation error of the bridge", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "event", seq: request.seq, json: JSON.stringify({ type: "log", data: {} }) });
    });

    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).followRun("r1"),
    ).rejects.toBeInstanceOf(ZodError);
  });
});

describe("desktop_bridge files", () => {
  test("readRunFile joins the chunks the back end streams", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({ kind: "chunk", seq: request.seq, bytes: new Uint8Array([40, 65]).buffer });
      reply({ kind: "chunk", seq: request.seq, bytes: new Uint8Array([41]).buffer });
      reply({ kind: "end", seq: request.seq });
    });

    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).readRunFile("r1", "clock.nwk"),
    ).resolves.toStrictEqual(new Uint8Array([40, 65, 41]));
    expect(backend.requests).toStrictEqual([{ kind: "read-file", seq: 0, id: "r1", path: "clock.nwk" }]);
  });

  test("readRunFile rejects with the error of the back end", async () => {
    const backend = fakeBackend((request, reply) => {
      reply({
        kind: "error",
        seq: request.seq,
        error: JSON.stringify({ code: "invalid_request", message: "file path `../x` must name a file", causes: [] }),
      });
    });

    await expect(
      createDesktopBridge(backend.connection, fakeShell({ picked: [] })).readRunFile("r1", "../x"),
    ).rejects.toMatchObject({
      response: { code: "invalid_request" },
    });
  });

  test("saveRunFile asks the shell to save the file and reports whether it was saved", async () => {
    const backend = fakeBackend(() => undefined);
    const shell = fakeShell({ picked: [], saved: { saved: true } });

    await expect(createDesktopBridge(backend.connection, shell).saveRunFile("r1", "out/a.nwk", "a.nwk")).resolves.toBe(
      true,
    );
    expect(shell.saves).toStrictEqual([{ id: "r1", path: "out/a.nwk", name: "a.nwk" }]);
  });

  test("saveRunArchive resolves false when the user cancels the save dialog", async () => {
    const backend = fakeBackend(() => undefined);
    const shell = fakeShell({ picked: [], saved: { saved: false } });

    await expect(createDesktopBridge(backend.connection, shell).saveRunArchive("r1", "run.zip")).resolves.toBe(false);
    expect(shell.saves).toStrictEqual([{ id: "r1", name: "run.zip" }]);
  });

  test("a save the back end refuses rejects with its typed error", async () => {
    const response = { code: "invalid_request", message: "file path `../x` must name a file", causes: [] };
    const backend = fakeBackend(() => undefined);
    const shell = fakeShell({ picked: [], saved: { error: JSON.stringify(response) } });

    await expect(createDesktopBridge(backend.connection, shell).saveRunFile("r1", "../x", "x")).rejects.toMatchObject({
      response,
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

describe("desktop_bridge window connection", () => {
  test("a port the preload posts to the window connects the back end", () => {
    const target = fakeWindow();
    const shell = fakeShell({ picked: [] });
    const connection = windowBackendConnection(target, shell);
    const endpoints: ClientEndpoint[] = [];

    connection.onEndpoint((endpoint) => {
      endpoints.push(endpoint);
    });
    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [fakePort()] });

    expect([shell.connections, endpoints.length]).toStrictEqual([1, 1]);
  });

  test("messages from another source or without a port are ignored", () => {
    const target = fakeWindow();
    const connection = windowBackendConnection(target, fakeShell({ picked: [] }));
    const endpoints: ClientEndpoint[] = [];

    connection.onEndpoint((endpoint) => {
      endpoints.push(endpoint);
    });
    target.emit({ source: {}, data: { channel: BACKEND_PORT_CHANNEL }, ports: [fakePort()] });
    target.emit({ source: target, data: { channel: BACKEND_PORT_CHANNEL }, ports: [] });
    target.emit({ source: target, data: "treetime", ports: [fakePort()] });

    expect(endpoints).toStrictEqual([]);
  });

  test("the stop message of the preload carries the reason", () => {
    const target = fakeWindow();
    const connection = windowBackendConnection(target, fakeShell({ picked: [] }));
    const reasons: string[] = [];

    connection.onStopped((reason) => {
      reasons.push(reason);
    });
    target.emit({ source: target, data: { channel: BACKEND_STOPPED_CHANNEL, reason: "crashed" }, ports: [] });

    expect(reasons).toStrictEqual(["crashed"]);
  });
});

type Responder = (request: BackendRequest, reply: (message: BackendReply) => void) => void;

interface FakeBackend {
  connection: BackendConnection;
  requests: BackendRequest[];
  calls(): unknown[];
  connect(): void;
  stop(reason: string): void;
}

function fakeBackend(respond: Responder, { connected = true } = {}): FakeBackend {
  const endpointListeners: Array<(endpoint: ClientEndpoint) => void> = [];
  const stoppedListeners: Array<(reason: string) => void> = [];
  let deliver: (reply: BackendReply) => void = () => undefined;

  const endpoint: ClientEndpoint = {
    post(request) {
      backend.requests.push(request);
      respond(request, (reply) => {
        deliver(reply);
      });
    },
    listen(listener) {
      deliver = listener;
    },
    onClose() {
      return undefined;
    },
  };

  const backend: FakeBackend = {
    requests: [],
    connection: {
      onEndpoint(listener) {
        endpointListeners.push(listener);

        if (connected) {
          listener(endpoint);
        }
      },
      onStopped(listener) {
        stoppedListeners.push(listener);
      },
    },
    calls: () => backend.requests.flatMap((request) => (request.kind === "call" ? [parseJson(request.request)] : [])),
    connect() {
      endpointListeners.forEach((listener) => {
        listener(endpoint);
      });
    },
    stop(reason) {
      stoppedListeners.forEach((listener) => {
        listener(reason);
      });
    },
  };

  return backend;
}

interface FakeShell extends DesktopShell {
  connections: number;
  saves: unknown[];
}

function fakeShell({ picked, saved = null }: { picked: unknown; saved?: unknown }): FakeShell {
  const shell: FakeShell = {
    connections: 0,
    saves: [],
    connectBackend() {
      shell.connections += 1;
    },
    pickFiles: () => Promise.resolve(picked),
    pathForFile: () => "",
    saveRunFile: (request: unknown) => {
      shell.saves.push(request);

      return Promise.resolve(saved);
    },
    saveRunArchive: (request: unknown) => {
      shell.saves.push(request);

      return Promise.resolve(saved);
    },
  };

  return shell;
}

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

function fakePort(): PortLike {
  return { postMessage: () => undefined, addEventListener: () => undefined, start: () => undefined };
}

function runEvent(seq: number, type: string, data: unknown): string {
  return JSON.stringify({ seq, time: "t", type, data });
}

function parseJson(text: string): unknown {
  const value: unknown = JSON.parse(text);

  return value;
}
