import { ApiError } from "@neherlab/app-contracts/client";
import type { PortMessage, PortReply } from "@neherlab/app-napi";
import { describe, expect, test, vi } from "vitest";

import { HostConnection, mainFetchPort, type MainPort } from "../host-port";
import { createPortFetch } from "../port-fetch";

const ORIGIN = "http://treetime.host";

describe("host_port adapter", () => {
  test("a main-process port carries requests out and replies in", () => {
    const port = fakeMainPort();
    const adapted = mainFetchPort(port);
    const received: PortReply[] = [];

    adapted.addEventListener("message", (event) => {
      received.push(event.data);
    });
    adapted.start();
    // oxlint-disable-next-line unicorn/require-post-message-target-origin -- a message port takes no target origin
    adapted.postMessage({ kind: "abort", seq: 1 });
    port.deliver({ kind: "end", seq: 1 });

    expect([port.started, port.posted, received]).toStrictEqual([
      true,
      [{ kind: "abort", seq: 1 }],
      [{ kind: "end", seq: 1 }],
    ]);
  });
});

describe("host_port connection", () => {
  test("a save during a restart waits for the port of the new back end", async () => {
    const connection = new HostConnection();
    const first = fakeMainPort();
    connection.connect(mainFetchPort(first));
    const portFetch = createPortFetch(connection);

    connection.stop({ reason: "the back end restarts", restarts: true });
    const saved = portFetch(`${ORIGIN}/api/runs/r1/save`, { method: "POST", body: "{}" });
    const second = fakeMainPort();
    connection.connect(mainFetchPort(second));
    await vi.waitFor(() => {
      expect(second.posted).toHaveLength(1);
    });
    second.deliver({ kind: "head", seq: 0, status: 200, headers: [] });
    second.deliver({ kind: "end", seq: 0 });

    await expect(saved.then((response) => response.status)).resolves.toBe(200);
    expect(first.posted).toStrictEqual([]);
  });

  test("a save after a failed start rejects at once with the start error", async () => {
    const connection = new HostConnection();
    connection.connect(mainFetchPort(fakeMainPort()));
    const portFetch = createPortFetch(connection);
    const error = { code: "internal_error", message: "When opening the runs folder '/r'", causes: [] } as const;

    connection.stop({ reason: "the back end could not start", restarts: false, error: { ...error, causes: [] } });

    await expect(portFetch(`${ORIGIN}/api/runs/r1/save`, { method: "POST", body: "{}" })).rejects.toStrictEqual(
      new ApiError(500, { ...error, causes: [] }),
    );
  });

  test("a listener added after the connection gets the current port", () => {
    const connection = new HostConnection();
    const port = mainFetchPort(fakeMainPort());
    const received: unknown[] = [];

    connection.connect(port);
    connection.onPort((next) => {
      received.push(next);
    });

    expect(received).toStrictEqual([port]);
  });
});

interface FakeMainPort extends MainPort {
  posted: PortMessage[];
  started: boolean;
  deliver(reply: PortReply): void;
}

function fakeMainPort(): FakeMainPort {
  const listeners: Array<(event: { data: PortReply }) => void> = [];

  const port: FakeMainPort = {
    posted: [],
    started: false,
    postMessage(message) {
      port.posted.push(message);
    },
    on(_event, listener) {
      listeners.push(listener);
    },
    start() {
      port.started = true;
    },
    deliver(reply) {
      listeners.forEach((listener) => {
        listener({ data: reply });
      });
    },
  };

  return port;
}
