import type { BackendStopped } from "@neherlab/app-ui/host";
import { describe, expect, test } from "vitest";

import { emit, handle, listenForPortRequests, sendPort, type IpcMainLike, type SenderEvent } from "../ipc-main";

const APP_URL = "app://treetime/";

const TRUSTED = { url: "app://treetime/index.html", parent: null };

const FOREIGN = { url: "https://example.org/", parent: null };

describe("ipc_main handle", () => {
  test("a trusted request reaches the handler parsed, and the reply goes back", async () => {
    const ipc = fakeIpcMain();
    const requests: unknown[] = [];

    handle(ipc, APP_URL, "pick-folder", (_sender, request) => {
      requests.push(request);

      return "/data/runs";
    });

    await expect(ipc.invoke("treetime:pick-folder", TRUSTED, { title: "Runs" })).resolves.toBe("/data/runs");
    expect(requests).toStrictEqual([{ title: "Runs" }]);
  });

  test("a request from a frame outside the application is refused before the handler runs", async () => {
    const ipc = fakeIpcMain();
    let called = false;

    handle(ipc, APP_URL, "restart-backend", () => {
      called = true;

      return undefined;
    });

    await expect(ipc.invoke("treetime:restart-backend", FOREIGN, undefined)).rejects.toThrow(
      "treetime:restart-backend refused a message from a frame outside the application",
    );
    expect(called).toBe(false);
  });

  test("a request that fails the schema of its channel is refused before the handler runs", async () => {
    const ipc = fakeIpcMain();
    let called = false;

    handle(ipc, APP_URL, "save-run", () => {
      called = true;

      return { kind: "canceled" as const };
    });

    await expect(
      ipc.invoke("treetime:save-run", TRUSTED, { id: "r1", name: "r1.zip", destination: "/etc/passwd" }),
    ).rejects.toMatchObject({ issues: [{ code: "unrecognized_keys", keys: ["destination"] }] });
    expect(called).toBe(false);
  });

  test("a theme outside the generated theme values is refused", async () => {
    const ipc = fakeIpcMain();

    handle(ipc, APP_URL, "set-native-theme", () => undefined);

    await expect(ipc.invoke("treetime:set-native-theme", TRUSTED, "purple")).rejects.toMatchObject({
      issues: [{ code: "invalid_union" }],
    });
  });
});

describe("ipc_main events and ports", () => {
  test("a port request from the application reaches the listener with its sender", () => {
    const ipc = fakeIpcMain();
    const senders: string[] = [];

    listenForPortRequests(ipc, APP_URL, (sender) => {
      senders.push(sender);
    });
    ipc.send("treetime:backend-port-request", TRUSTED);
    ipc.send("treetime:backend-port-request", FOREIGN);

    expect(senders).toStrictEqual(["window"]);
  });

  test("an event goes to the prefixed channel after its schema check", () => {
    const contents = fakeContents();
    const stop: BackendStopped = { reason: "crashed", restarts: false };

    emit(contents, "backend-stopped", stop);

    expect(contents.sent).toStrictEqual([["treetime:backend-stopped", stop]]);
  });

  test("a port goes to the window on the port channel", () => {
    const contents = fakeContents();

    sendPort(contents, "port");

    expect(contents.posted).toStrictEqual([["treetime:backend-port", null, ["port"]]]);
  });
});

interface FakeIpcMain extends IpcMainLike<string> {
  invoke(channel: string, frame: SenderEvent<string>["senderFrame"], request: unknown): Promise<unknown>;
  send(channel: string, frame: SenderEvent<string>["senderFrame"]): void;
}

function fakeIpcMain(): FakeIpcMain {
  const handlers = new Map<string, (event: SenderEvent<string>, ...args: unknown[]) => unknown>();
  const listeners = new Map<string, (event: SenderEvent<string>, ...args: unknown[]) => void>();

  return {
    handle(channel, listener) {
      handlers.set(channel, listener);
    },
    on(channel, listener) {
      listeners.set(channel, listener);
    },
    invoke(channel, senderFrame, request) {
      return Promise.resolve(handlers.get(channel)?.({ sender: "window", senderFrame }, request));
    },
    send(channel, senderFrame) {
      listeners.get(channel)?.({ sender: "window", senderFrame });
    },
  };
}

function fakeContents() {
  const sent: unknown[][] = [];
  const posted: unknown[][] = [];

  return {
    sent,
    posted,
    send(channel: string, ...args: unknown[]) {
      sent.push([channel, ...args]);
    },
    postMessage(channel: string, message: null, transfer: string[]) {
      posted.push([channel, message, transfer]);
    },
  };
}
