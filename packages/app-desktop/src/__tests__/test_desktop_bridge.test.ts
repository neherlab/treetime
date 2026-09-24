import { CancelledError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { createDesktopBridge, type IpcRendererLike } from "../desktop-bridge";

type Listener = (event: unknown, ...args: unknown[]) => void;

type OnInvoke = (ipc: FakeIpc, channel: string, argsJson: unknown) => Promise<unknown>;

interface FakeIpc extends IpcRendererLike {
  emit(channel: string, data: unknown): void;
  handlerCount(channel: string): number;
  readonly sent: string[];
}

function makeFakeIpc(onInvoke: OnInvoke): FakeIpc {
  const handlers = new Map<string, Listener[]>();
  const sent: string[] = [];

  const fake: FakeIpc = {
    invoke: (channel, ...args) => onInvoke(fake, channel, args[0]),
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
    send(channel) {
      sent.push(channel);
    },
    emit(channel, data) {
      for (const listener of handlers.get(channel) ?? []) {
        listener(undefined, data);
      }
    },
    handlerCount(channel) {
      return (handlers.get(channel) ?? []).length;
    },
    sent,
  };

  return fake;
}

describe("desktop_bridge query path", () => {
  test("version parses an object result", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve({ version: "1.0.0" })));
    await expect(bridge.version()).resolves.toStrictEqual({ version: "1.0.0" });
  });

  test("version parses a JSON string result", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.resolve(JSON.stringify({ version: "2.0.0" }))));
    await expect(bridge.version()).resolves.toStrictEqual({ version: "2.0.0" });
  });
});

describe("desktop_bridge streaming command path", () => {
  test("progress and log events reach the caller and the result validates", async () => {
    const fake = makeFakeIpc((ipc, channel) => {
      if (channel === "treetime:ancestral") {
        ipc.emit("treetime:progress", { stage: "infer", fraction: 1, message: "done" });
        ipc.emit("treetime:log", { level: "Info", message: "ok" });

        return Promise.resolve({ model_name: "JC69" });
      }

      return Promise.resolve(undefined);
    });

    const received: string[] = [];
    const bridge = createDesktopBridge(fake);

    const result = await bridge.ancestral(
      { tree: "t", outdir: "o" },
      {
        onProgress: (e) => {
          received.push(e.stage);
        },
      },
    );

    expect(result).toStrictEqual({ model_name: "JC69" });
    expect(received).toStrictEqual(["infer"]);
    expect(fake.handlerCount("treetime:progress")).toBe(0);
    expect(fake.handlerCount("treetime:log")).toBe(0);
  });

  test("an IPC error mentioning cancellation surfaces as CancelledError", async () => {
    const bridge = createDesktopBridge(makeFakeIpc(() => Promise.reject(new Error("Computation cancelled by client"))));
    await expect(bridge.clock({ dates: "d", outdir: "o" })).rejects.toBeInstanceOf(CancelledError);
  });

  test("aborting the signal sends the cancel channel and cleans up listeners", async () => {
    let resolveInvoke: (value: unknown) => void = () => undefined;

    const fake = makeFakeIpc(
      () =>
        new Promise<unknown>((resolve) => {
          resolveInvoke = resolve;
        }),
    );

    const controller = new AbortController();
    const bridge = createDesktopBridge(fake);
    const pending = bridge.ancestral({ tree: "t", outdir: "o" }, { signal: controller.signal });

    controller.abort();
    expect(fake.sent).toContain("treetime:cancel");

    resolveInvoke({ model_name: "JC69" });
    await expect(pending).resolves.toStrictEqual({ model_name: "JC69" });
    expect(fake.handlerCount("treetime:progress")).toBe(0);
  });
});
