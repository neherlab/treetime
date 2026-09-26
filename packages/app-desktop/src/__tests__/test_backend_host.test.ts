import { BridgeError } from "@neherlab/app-contracts";
import type { PortMessage, PortReply, PortRequest, Subscription } from "@neherlab/app-napi";
import { describe, expect, test } from "vitest";

import { saveRunFiles, serveBackend, serveFetch, type AddonBackend } from "../backend-host";
import {
  portEndpoint,
  zBackendReply,
  zBackendRequest,
  type BackendReply,
  type BackendRequest,
  type FetchEndpoint,
  type HostEndpoint,
  type PortLike,
} from "../backend-protocol";
import { createDesktopBridge, type BackendConnection } from "../desktop-bridge";

type EventCallback = (err: Error | null, eventJson: string) => void;

describe("backend_host requests", () => {
  test("a call answers with the JSON result of the addon", async () => {
    const host = hostWith({ call: (request) => Promise.resolve(`{"echo":${request}}`) });

    host.send({ kind: "call", seq: 4, request: '{"operation":"version","args":{}}' });

    await expect(host.replies(1)).resolves.toStrictEqual([
      { kind: "result", seq: 4, json: '{"echo":{"operation":"version","args":{}}}' },
    ]);
  });

  test("a rejected call answers with the error message of the addon", async () => {
    const response = '{"code":"conflict","message":"run `a` has already started","causes":[]}';
    const host = hostWith({ call: () => Promise.reject(new Error(response)) });

    host.send({ kind: "call", seq: 1, request: "{}" });

    await expect(host.replies(1)).resolves.toStrictEqual([{ kind: "error", seq: 1, error: response }]);
  });

  test("a subscription forwards events until it is unsubscribed", async () => {
    let forward: EventCallback = () => undefined;
    const unsubscribed: number[] = [];

    const host = hostWith({
      subscribe: (_id, _from, onEvent) => {
        forward = onEvent;

        return subscription(() => {
          unsubscribed.push(1);
        });
      },
    });

    host.send({ kind: "subscribe", seq: 2, id: "r1", from: 0 });
    forward(null, '{"seq":0}');
    host.send({ kind: "unsubscribe", seq: 2 });

    await expect(host.replies(1)).resolves.toStrictEqual([{ kind: "event", seq: 2, json: '{"seq":0}' }]);
    expect(unsubscribed).toStrictEqual([1]);
  });

  test("a subscription the addon refuses answers with its error", async () => {
    const host = hostWith({
      subscribe: () => {
        throw new Error('{"code":"not_found","message":"no run `r9`","causes":[]}');
      },
    });

    host.send({ kind: "subscribe", seq: 3, id: "r9", from: 0 });

    await expect(host.replies(1)).resolves.toStrictEqual([
      { kind: "error", seq: 3, error: '{"code":"not_found","message":"no run `r9`","causes":[]}' },
    ]);
  });

  test("closing the port ends every subscription of the renderer", () => {
    const unsubscribed: string[] = [];

    const host = hostWith({
      subscribe: (id) =>
        subscription(() => {
          unsubscribed.push(id);
        }),
    });

    host.send({ kind: "subscribe", seq: 1, id: "r1", from: 0 });
    host.send({ kind: "subscribe", seq: 2, id: "r2", from: 0 });
    host.close();

    expect(unsubscribed).toStrictEqual(["r1", "r2"]);
  });
});

describe("backend_host fetch relay", () => {
  test("a request goes to the addon and every reply of the addon goes back to the port", () => {
    const received: PortRequest[] = [];

    const host = fetchHostWith((request, onReply) => {
      received.push(request);
      onReply({
        kind: "head",
        seq: request.seq,
        status: 200,
        headers: [{ name: "content-type", value: "text/plain" }],
      });
      onReply({ kind: "chunk", seq: request.seq, data: new Uint8Array([104, 105]) });
      onReply({ kind: "end", seq: request.seq });

      return { abort: () => undefined };
    });

    host.send({ kind: "request", request: REQUEST });

    expect([received, host.replies]).toStrictEqual([
      [REQUEST],
      [
        { kind: "head", seq: 3, status: 200, headers: [{ name: "content-type", value: "text/plain" }] },
        { kind: "chunk", seq: 3, data: new Uint8Array([104, 105]) },
        { kind: "end", seq: 3 },
      ],
    ]);
  });

  test("an abort message aborts the open exchange of its sequence number only", () => {
    const aborted: number[] = [];

    const host = fetchHostWith((request) => ({
      abort: () => {
        aborted.push(request.seq);
      },
    }));

    host.send({ kind: "request", request: { ...REQUEST, seq: 1 } });
    host.send({ kind: "request", request: { ...REQUEST, seq: 2 } });
    host.send({ kind: "abort", seq: 2 });
    host.send({ kind: "abort", seq: 2 });

    expect(aborted).toStrictEqual([2]);
  });

  test("an ended exchange is not aborted by a late abort message or by closing the port", () => {
    const aborted: number[] = [];

    const host = fetchHostWith((request, onReply) => {
      onReply({ kind: "end", seq: request.seq });

      return {
        abort: () => {
          aborted.push(request.seq);
        },
      };
    });

    host.send({ kind: "request", request: REQUEST });
    host.send({ kind: "abort", seq: REQUEST.seq });
    host.close();

    expect(aborted).toStrictEqual([]);
  });

  test("closing the port aborts every open exchange", () => {
    const aborted: number[] = [];

    const host = fetchHostWith((request) => ({
      abort: () => {
        aborted.push(request.seq);
      },
    }));

    host.send({ kind: "request", request: { ...REQUEST, seq: 1 } });
    host.send({ kind: "request", request: { ...REQUEST, seq: 2 } });
    host.close();

    expect(aborted).toStrictEqual([1, 2]);
  });

  test("a request the addon refuses answers with an invalid request error", () => {
    const host = fetchHostWith(() => {
      throw new Error("Failed to convert JavaScript value `Undefined` into rust type `String` on PortRequest.url");
    });

    host.send({ kind: "request", request: REQUEST });

    expect(host.replies).toStrictEqual([
      {
        kind: "error",
        seq: 3,
        error: {
          code: "invalid_request",
          message: "Failed to convert JavaScript value `Undefined` into rust type `String` on PortRequest.url",
          causes: [],
        },
      },
    ]);
  });
});

describe("backend_host saves", () => {
  test("a file save writes through the addon and reports it", async () => {
    const saves: string[][] = [];

    const reply = await saveRunFiles(
      addonWith({
        saveRunFile: ({ id, path, destination }) => {
          saves.push([id, path, destination]);

          return Promise.resolve();
        },
      }),
      { kind: "save-file", seq: 3, request: { id: "r1", path: "a.nwk", destination: "/home/user/a.nwk" } },
    );

    expect([reply, saves]).toStrictEqual([{ kind: "saved", seq: 3 }, [["r1", "a.nwk", "/home/user/a.nwk"]]]);
  });

  test("an archive save that fails reports the error of the addon", async () => {
    const reply = await saveRunFiles(
      addonWith({ saveRunArchive: () => Promise.reject(new Error("When saving '/ro/r1.zip': permission denied")) }),
      { kind: "save-archive", seq: 4, request: { id: "r1", destination: "/ro/r1.zip" } },
    );

    expect(reply).toStrictEqual({ kind: "error", seq: 4, error: "When saving '/ro/r1.zip': permission denied" });
  });
});

describe("backend_host with the renderer bridge over a message channel", () => {
  test("the bridge reaches the addon and receives its typed errors", async () => {
    const channel = new MessageChannel();

    const addon: AddonBackend = {
      call: (request) =>
        request.includes('"version"')
          ? Promise.resolve('{"version":"2.0.0"}')
          : Promise.reject(new Error('{"code":"not_found","message":"no run `r9`","causes":[]}')),
      subscribe: () => subscription(() => undefined),
      saveRunFile: () => Promise.resolve(),
      saveRunArchive: () => Promise.resolve(),
      fetch: () => ({ abort: () => undefined }),
    };

    const hostPort = nodePort(channel.port1);
    const rendererPort = nodePort(channel.port2);
    serveBackend(portEndpoint<BackendRequest, BackendReply>(hostPort, zBackendRequest), addon);

    const bridge = createDesktopBridge(connected(rendererPort), {
      connectBackend: () => undefined,
      pickFiles: () => Promise.resolve([]),
      pathForFile: () => "",
      saveRunFile: () => Promise.resolve({ saved: true }),
      saveRunArchive: () => Promise.resolve({ saved: true }),
    });

    try {
      await expect(bridge.version()).resolves.toStrictEqual({ version: "2.0.0" });
      await expect(bridge.getRun("r9")).rejects.toBeInstanceOf(BridgeError);
    } finally {
      channel.port1.close();
      channel.port2.close();
    }
  });
});

const REQUEST: PortRequest = { seq: 3, method: "GET", url: "/api/version", headers: [] };

function fetchHostWith(fetch: AddonBackend["fetch"]) {
  const listeners: Array<(message: PortMessage) => void> = [];
  const closeListeners: Array<() => void> = [];
  const replies: PortReply[] = [];

  const endpoint: FetchEndpoint = {
    post(reply) {
      replies.push(reply);
    },
    listen(listener) {
      listeners.push(listener);
    },
    onClose(listener) {
      closeListeners.push(listener);
    },
  };

  serveFetch(endpoint, { fetch });

  return {
    replies,
    send(message: PortMessage) {
      listeners.forEach((listener) => {
        listener(message);
      });
    },
    close() {
      closeListeners.forEach((listener) => {
        listener();
      });
    },
  };
}

function hostWith(overrides: Partial<AddonBackend>) {
  const listeners: Array<(request: BackendRequest) => void> = [];
  const closeListeners: Array<() => void> = [];
  const received: BackendReply[] = [];
  const waiters: Array<() => void> = [];

  const endpoint: HostEndpoint = {
    post(message) {
      received.push(message);
      waiters.splice(0).forEach((wake) => {
        wake();
      });
    },
    listen(listener) {
      listeners.push(listener);
    },
    onClose(listener) {
      closeListeners.push(listener);
    },
  };

  serveBackend(endpoint, addonWith(overrides));

  return {
    send(request: BackendRequest) {
      listeners.forEach((listener) => {
        listener(request);
      });
    },
    close() {
      closeListeners.forEach((listener) => {
        listener();
      });
    },
    async replies(count: number): Promise<BackendReply[]> {
      while (received.length < count) {
        await new Promise<void>((wake) => {
          waiters.push(wake);
        });
      }

      return received;
    },
  };
}

function addonWith(overrides: Partial<AddonBackend>): AddonBackend {
  return {
    call: () => Promise.reject(new Error("call is not part of this test")),
    subscribe: () => {
      throw new Error("subscribe is not part of this test");
    },
    saveRunFile: () => Promise.reject(new Error("saveRunFile is not part of this test")),
    saveRunArchive: () => Promise.reject(new Error("saveRunArchive is not part of this test")),
    fetch: () => {
      throw new Error("fetch is not part of this test");
    },
    ...overrides,
  };
}

function subscription(unsubscribe: () => void): Pick<Subscription, "unsubscribe"> {
  return { unsubscribe };
}

function connected(port: PortLike): BackendConnection {
  return {
    onEndpoint: (listener) => {
      listener(portEndpoint<BackendReply, BackendRequest>(port, zBackendReply));
    },
    onStopped: () => undefined,
  };
}

function nodePort(port: InstanceType<typeof MessageChannel>["port1"]): PortLike {
  return {
    postMessage: (message) => {
      port.postMessage(message);
    },
    addEventListener: (type, listener) => {
      port.on(type, (data) => {
        listener({ data });
      });
    },
    start: () => {
      port.start();
    },
  };
}
