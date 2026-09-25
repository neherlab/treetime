import { BridgeError, CancelledError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { createDesktopBridge, createLocalFiles, RUN_EVENT_CHANNEL, type IpcRendererLike } from "../desktop-bridge";

type Listener = (event: unknown, ...args: unknown[]) => void;

type OnInvoke = (ipc: FakeIpc, channel: string, args: unknown[]) => Promise<unknown>;

interface FakeIpc extends IpcRendererLike {
  emit(channel: string, ...args: unknown[]): void;
  handlerCount(channel: string): number;
  sent: unknown[][];
}

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

function makeFakeIpc(onInvoke: OnInvoke): FakeIpc {
  const handlers = new Map<string, Listener[]>();

  const fake: FakeIpc = {
    sent: [],
    invoke: (channel, ...args) =>
      onInvoke(fake, channel, args).then(
        (value) => ({ ok: true, value }),
        (error: unknown) => ({ ok: false, error: error instanceof Error ? error.message : String(error) }),
      ),
    on(channel, listener) {
      handlers.set(channel, [...(handlers.get(channel) ?? []), listener]);
    },
    removeListener(channel, listener) {
      handlers.set(
        channel,
        (handlers.get(channel) ?? []).filter((l) => l !== listener),
      );
    },
    send(channel, ...args) {
      fake.sent.push([channel, ...args]);
    },
    emit(channel, ...args) {
      for (const listener of handlers.get(channel) ?? []) {
        listener(undefined, ...args);
      }
    },
    handlerCount(channel) {
      return (handlers.get(channel) ?? []).length;
    },
  };

  return fake;
}

function runEvent(seq: number, type: string, data: unknown): string {
  return JSON.stringify({ seq, time: "t", type, data });
}

describe("desktop_bridge queries and requests", () => {
  test("version parses a JSON string result", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve(JSON.stringify({ version: "2.0.0" }))));
    await expect(bridge.version()).resolves.toStrictEqual({ version: "2.0.0" });
  });

  test("checkConfig sends the check-config operation with its request", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(
          JSON.stringify({
            status: "valid",
            command: "clock",
            config: {},
            code: { command_line: [], command_line_text: "", yaml: [], yaml_text: "" },
            checks: [],
          }),
        );
      }),
    );

    await bridge.checkConfig({ command: "clock", text: "tree: t" });
    expect(captured).toStrictEqual([
      "treetime:call",
      JSON.stringify({ operation: "check-config", args: { request: { command: "clock", text: "tree: t" } } }),
    ]);
  });

  test("runConfig sends the run-config operation with its request", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(
          JSON.stringify({
            status: "valid",
            command: "clock",
            config: {},
            code: { command_line: [], command_line_text: "", yaml: [], yaml_text: "" },
            checks: [],
          }),
        );
      }),
    );

    await bridge.runConfig({ command: "clock", config: { tree: "t" } });
    expect(captured).toStrictEqual([
      "treetime:call",
      JSON.stringify({ operation: "run-config", args: { request: { command: "clock", config: { tree: "t" } } } }),
    ]);
  });

  test("startRun sends the start-run operation with the replacement configuration", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(JSON.stringify(RECORD));
      }),
    );

    await bridge.startRun("r1", { config: { tree: "/data/t.nwk" } });
    expect(captured).toStrictEqual([
      "treetime:call",
      JSON.stringify({ operation: "start-run", args: { id: "r1", request: { config: { tree: "/data/t.nwk" } } } }),
    ]);
  });

  test("readRunFile returns the bytes the main process sends", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve(new Uint8Array([40, 65, 41]))));
    await expect(bridge.readRunFile("r1", "clock.nwk")).resolves.toStrictEqual(new Uint8Array([40, 65, 41]));
  });

  test("readRunFile rejects a result without bytes", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve("not bytes")));
    await expect(bridge.readRunFile("r1", "clock.nwk")).rejects.toBeInstanceOf(TypeError);
  });

  test("uploadInput rejects because the desktop app reads local paths", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve(null)));
    await expect(bridge.uploadInput("r1", "t.nwk", new Blob(["x"]))).rejects.toThrow(
      "the desktop application reads inputs from local file paths; name the files in the run configuration",
    );
  });
});

describe("desktop_bridge errors", () => {
  test("a typed error of the back end rejects with a BridgeError carrying its class and causes", async () => {
    const response = { code: "not_found", message: "When reading run `r9`", causes: ["no run with id `r9`"] };
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.reject(new Error(JSON.stringify(response)))));

    const error = await bridge.getRun("r9").catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(BridgeError);
    expect(error).toMatchObject({ message: "When reading run `r9`: no run with id `r9`", response });
  });

  test("an untyped error rejects as an internal error with its message", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.reject(new Error("addon failed to load"))));

    await expect(bridge.version()).rejects.toMatchObject({
      response: { code: "internal_error", message: "addon failed to load", causes: [] },
    });
  });

  test("a reply without the result envelope rejects", async () => {
    const ipc: IpcRendererLike = {
      invoke: () => Promise.resolve("{}"),
      on: () => undefined,
      removeListener: () => undefined,
      send: () => undefined,
    };

    await expect(createDesktopBridge(ipc).version()).rejects.toThrow("the main process sent a reply without a result");
  });
});

describe("desktop_bridge run events", () => {
  test("followRun subscribes and resolves at the terminal event of its own subscription", async () => {
    let subscription: unknown[] = [];

    const ipc = makeFakeIpc((fake, channel, args) => {
      if (channel === "treetime:runs:subscribe") {
        subscription = args;
        fake.emit(RUN_EVENT_CHANNEL, "other", runEvent(3, "terminal", { status: "cancelled", job_id: "x" }));
        fake.emit(RUN_EVENT_CHANNEL, "sub-1", runEvent(3, "log", { level: "info", message: "hello" }));
        fake.emit(RUN_EVENT_CHANNEL, "sub-1", runEvent(4, "terminal", { status: "cancelled", job_id: "r1" }));
      }

      return Promise.resolve(undefined);
    });

    const bridge = createDesktopBridge(ipc, () => "sub-1");
    const seen: number[] = [];

    const terminal = await bridge.followRun("r1", {
      from: 3,
      onEvent: (event) => {
        seen.push(event.seq);
      },
    });

    expect(subscription).toStrictEqual(["sub-1", "r1", 3]);
    expect(terminal).toStrictEqual({ status: "cancelled", job_id: "r1" });
    expect(seen).toStrictEqual([3, 4]);
    expect(ipc.handlerCount(RUN_EVENT_CHANNEL)).toBe(0);
    expect(ipc.sent).toStrictEqual([["treetime:runs:unsubscribe", "sub-1"]]);
  });

  test("aborting followRun unsubscribes the run events in the main process", async () => {
    const ipc = makeFakeIpc(() => Promise.resolve(undefined));
    const controller = new AbortController();
    const bridge = createDesktopBridge(ipc, () => "sub-1");

    const following = bridge.followRun("r1", { signal: controller.signal });
    await Promise.resolve();
    controller.abort();

    await expect(following).rejects.toThrow("the event stream of run r1 ended without a terminal event");
    expect(ipc.handlerCount(RUN_EVENT_CHANNEL)).toBe(0);
    expect(ipc.sent).toStrictEqual([["treetime:runs:unsubscribe", "sub-1"]]);
  });

  test("a command creates a run and resolves with the outcome", async () => {
    const ipc = makeFakeIpc((fake, channel) => {
      if (channel === "treetime:call") {
        return Promise.resolve(JSON.stringify(RECORD));
      }

      if (channel === "treetime:runs:subscribe") {
        fake.emit(RUN_EVENT_CHANNEL, "sub-1", runEvent(0, "started", { job_id: "r1", command: "clock" }));
        fake.emit(RUN_EVENT_CHANNEL, "sub-1", runEvent(1, "terminal", { status: "ok", job_id: "r1", result: OUTCOME }));
      }

      return Promise.resolve(undefined);
    });

    const bridge = createDesktopBridge(ipc, () => "sub-1");
    await expect(bridge.clock({ tree: "t.nwk" })).resolves.toStrictEqual(OUTCOME);
  });

  test("a cancelled run rejects the command with a CancelledError", async () => {
    const ipc = makeFakeIpc((fake, channel) => {
      if (channel === "treetime:call") {
        return Promise.resolve(JSON.stringify(RECORD));
      }

      if (channel === "treetime:runs:subscribe") {
        fake.emit(RUN_EVENT_CHANNEL, "sub-1", runEvent(0, "terminal", { status: "cancelled", job_id: "r1" }));
      }

      return Promise.resolve(undefined);
    });

    const bridge = createDesktopBridge(ipc, () => "sub-1");
    await expect(bridge.clock({ tree: "t.nwk" })).rejects.toBeInstanceOf(CancelledError);
  });

  test("a failed subscription rejects and removes its listener", async () => {
    const ipc = makeFakeIpc((_fake, channel) =>
      channel === "treetime:runs:subscribe" ? Promise.reject(new Error("no run with id `r9`")) : Promise.resolve(null),
    );

    const bridge = createDesktopBridge(ipc, () => "sub-1");
    await expect(bridge.followRun("r9")).rejects.toThrow("no run with id `r9`");
    expect(ipc.handlerCount(RUN_EVENT_CHANNEL)).toBe(0);
  });
});

describe("desktop_bridge local files", () => {
  test("pickFiles asks the main process through the pick-files channel", async () => {
    let captured: unknown[] = [];

    const files = createLocalFiles(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(["/data/tree.nwk"]);
      }),
      () => "",
    );

    await expect(files.pickFiles({ title: "Tree", extensions: ["nwk"], multiple: false })).resolves.toStrictEqual([
      "/data/tree.nwk",
    ]);
    expect(captured).toStrictEqual([
      "treetime:pick-files",
      JSON.stringify({ title: "Tree", extensions: ["nwk"], multiple: false }),
    ]);
  });

  test("a malformed pick result rejects", async () => {
    const files = createLocalFiles(
      makeFakeIpc(() => Promise.resolve([1])),
      () => "",
    );

    await expect(files.pickFiles({ title: "Tree", extensions: [], multiple: false })).rejects.toBeInstanceOf(Error);
  });
});
