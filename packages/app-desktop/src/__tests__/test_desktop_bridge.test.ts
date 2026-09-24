import { CancelledError, CommandError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { createDesktopBridge, JOB_EVENT_CHANNEL, type IpcRendererLike } from "../desktop-bridge";

type Listener = (event: unknown, ...args: unknown[]) => void;

type OnInvoke = (ipc: FakeIpc, channel: string, args: unknown[]) => Promise<unknown>;

interface FakeIpc extends IpcRendererLike {
  emit(channel: string, ...args: unknown[]): void;
  handlerCount(channel: string): number;
  readonly sent: unknown[][];
}

const OUTCOME = { command: "ancestral", output_files: ["out/ancestral.nwk"] };

function makeFakeIpc(onInvoke: OnInvoke): FakeIpc {
  const handlers = new Map<string, Listener[]>();
  const sent: unknown[][] = [];

  const fake: FakeIpc = {
    invoke: (channel, ...args) => onInvoke(fake, channel, args),
    on(channel, listener) {
      const list = handlers.get(channel) ?? [];
      list.push(listener);
      handlers.set(channel, list);
    },
    removeListener(channel, listener) {
      handlers.set(
        channel,
        (handlers.get(channel) ?? []).filter((l) => l !== listener),
      );
    },
    send(channel, ...args) {
      sent.push([channel, ...args]);
    },
    emit(channel, ...args) {
      for (const listener of handlers.get(channel) ?? []) {
        listener(undefined, ...args);
      }
    },
    handlerCount(channel) {
      return (handlers.get(channel) ?? []).length;
    },
    sent,
  };

  return fake;
}

function jobIds(...ids: string[]): () => string {
  let next = 0;

  return () => ids[next++] ?? "unexpected";
}

describe("desktop_bridge query and request paths", () => {
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

    await expect(bridge.checkConfig({ command: "clock", text: "tree: t" })).resolves.toStrictEqual({
      status: "valid",
      config: {},
    });
    expect(captured).toStrictEqual(["treetime:check-config", JSON.stringify({ command: "clock", text: "tree: t" })]);
  });
});

describe("desktop_bridge command path", () => {
  test("events of this job reach the caller, events of other jobs do not", async () => {
    let captured: unknown[] = [];

    const fake = makeFakeIpc((ipc, channel, args) => {
      captured = [channel, ...args];
      ipc.emit(
        JOB_EVENT_CHANNEL,
        "other",
        JSON.stringify({ type: "progress", data: { stage: "x", fraction: 0, message: "" } }),
      );
      ipc.emit(
        JOB_EVENT_CHANNEL,
        "job-1",
        JSON.stringify({ type: "started", data: { job_id: "job-1", command: "ancestral" } }),
      );
      ipc.emit(
        JOB_EVENT_CHANNEL,
        "job-1",
        JSON.stringify({ type: "progress", data: { stage: "infer", fraction: 1, message: "" } }),
      );

      return Promise.resolve(JSON.stringify({ status: "ok", job_id: "job-1", result: OUTCOME }));
    });

    const received: string[] = [];
    const bridge = createDesktopBridge(fake, jobIds("job-1"));

    const result = await bridge.ancestral(
      { tree: "t.nwk" },
      {
        onStarted: (jobId) => {
          received.push(`started ${jobId}`);
        },
        onProgress: (event) => {
          received.push(event.stage);
        },
      },
    );

    expect(result).toStrictEqual(OUTCOME);
    expect(received).toStrictEqual(["started job-1", "infer"]);
    expect(captured).toStrictEqual(["treetime:run", "job-1", "ancestral", JSON.stringify({ tree: "t.nwk" })]);
    expect(fake.handlerCount(JOB_EVENT_CHANNEL)).toBe(0);
  });

  test("an error terminal event rejects with a CommandError", async () => {
    const terminal = { status: "error", job_id: "job-2", message: "bad", causes: [] };

    const bridge = createDesktopBridge(
      makeFakeIpc(() => Promise.resolve(JSON.stringify(terminal))),
      jobIds("job-2"),
    );

    await expect(bridge.clock({ tree: "t.nwk" })).rejects.toBeInstanceOf(CommandError);
  });

  test("aborting sends the job id on the cancel channel and the cancelled terminal event rejects", async () => {
    let resolveInvoke: (value: unknown) => void = () => undefined;

    const fake = makeFakeIpc(
      () =>
        new Promise<unknown>((resolve) => {
          resolveInvoke = resolve;
        }),
    );

    const controller = new AbortController();
    const bridge = createDesktopBridge(fake, jobIds("job-3"));
    const pending = bridge.ancestral({ tree: "t.nwk" }, { signal: controller.signal });

    controller.abort();
    expect(fake.sent).toStrictEqual([["treetime:cancel", "job-3"]]);

    resolveInvoke(JSON.stringify({ status: "cancelled", job_id: "job-3" }));
    await expect(pending).rejects.toBeInstanceOf(CancelledError);
    expect(fake.handlerCount(JOB_EVENT_CHANNEL)).toBe(0);
  });

  test("a signal aborted before the call rejects without starting a job", async () => {
    const invoked: string[] = [];

    const fake = makeFakeIpc((_ipc, channel) => {
      invoked.push(channel);

      return Promise.resolve(undefined);
    });

    const controller = new AbortController();
    controller.abort();
    const bridge = createDesktopBridge(fake, jobIds("job-4"));

    await expect(bridge.prune({ tree: "t.nwk" }, { signal: controller.signal })).rejects.toBeInstanceOf(CancelledError);
    expect(invoked).toStrictEqual([]);
  });
});
