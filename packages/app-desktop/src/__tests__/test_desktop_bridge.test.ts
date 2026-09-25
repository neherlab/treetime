import { CancelledError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { createDesktopBridge, LocalInputsError, RUN_EVENT_CHANNEL, type IpcRendererLike } from "../desktop-bridge";

type Listener = (event: unknown, ...args: unknown[]) => void;

type OnInvoke = (ipc: FakeIpc, channel: string, args: unknown[]) => Promise<unknown>;

interface FakeIpc extends IpcRendererLike {
  emit(channel: string, ...args: unknown[]): void;
  handlerCount(channel: string): number;
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
    invoke: (channel, ...args) => onInvoke(fake, channel, args),
    on(channel, listener) {
      handlers.set(channel, [...(handlers.get(channel) ?? []), listener]);
    },
    removeListener(channel, listener) {
      handlers.set(
        channel,
        (handlers.get(channel) ?? []).filter((l) => l !== listener),
      );
    },
    send() {
      return undefined;
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

  test("checkConfig sends the request as JSON to the check-config channel", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(JSON.stringify({ status: "valid", config: {} }));
      }),
    );

    await bridge.checkConfig({ command: "clock", text: "tree: t" });
    expect(captured).toStrictEqual(["treetime:check-config", JSON.stringify({ command: "clock", text: "tree: t" })]);
  });

  test("runConfig sends the request as JSON to the run-config channel", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(JSON.stringify({ status: "valid", config: {} }));
      }),
    );

    await bridge.runConfig({ command: "clock", config: { tree: "t" } });
    expect(captured).toStrictEqual(["treetime:run-config", JSON.stringify({ command: "clock", config: { tree: "t" } })]);
  });

  test("startRun sends the replacement configuration as JSON", async () => {
    let captured: unknown[] = [];

    const bridge = createDesktopBridge(
      makeFakeIpc((_ipc, channel, args) => {
        captured = [channel, ...args];

        return Promise.resolve(JSON.stringify(RECORD));
      }),
    );

    await bridge.startRun("r1", { config: { tree: "/data/t.nwk" } });
    expect(captured).toStrictEqual(["treetime:runs:start", "r1", JSON.stringify({ tree: "/data/t.nwk" })]);
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
    await expect(bridge.uploadInput("r1", "t.nwk", new Blob(["x"]))).rejects.toBeInstanceOf(LocalInputsError);
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
  });

  test("a command creates a run and resolves with the outcome", async () => {
    const ipc = makeFakeIpc((fake, channel) => {
      if (channel === "treetime:runs:create") {
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
      if (channel === "treetime:runs:create") {
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
