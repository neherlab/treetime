import type { PortMessage, PortReply, PortRequest } from "@neherlab/app-napi";
import { describe, expect, test } from "vitest";

import { serveFetch, type AddonBackend } from "../backend-host";
import type { FetchEndpoint } from "../backend-protocol";
import { createPortFetch, type FetchPort } from "../port-fetch";

describe("backend_host fetch relay", () => {
  test("a request goes to the addon and every reply of the addon goes back to the port", () => {
    const received: PortRequest[] = [];

    const host = fetchHostWith((request, _scope, onReply) => {
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

    const host = fetchHostWith((request, _scope, onReply) => {
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

  test("a request the addon refuses answers with the error response of the addon", () => {
    const rejected: Array<[number, string]> = [];

    const host = fetchHostWith(
      () => {
        throw new Error("Failed to convert JavaScript value `Undefined` into rust type `String` on PortRequest.url");
      },
      (seq, message) => {
        rejected.push([seq, message]);

        return [{ kind: "end", seq }];
      },
    );

    host.send({ kind: "request", request: REQUEST });

    expect([rejected, host.replies]).toStrictEqual([
      [[3, "Failed to convert JavaScript value `Undefined` into rust type `String` on PortRequest.url"]],
      [{ kind: "end", seq: 3 }],
    ]);
  });

  test("the scope of the port reaches the addon with every request", () => {
    const scopes: string[] = [];

    const host = fetchHostWith((_request, scope) => {
      scopes.push(scope);

      return { abort: () => undefined };
    });

    host.send({ kind: "request", request: REQUEST });

    expect(scopes).toStrictEqual(["renderer"]);
  });

  test("a reset reply ends the exchange, so a late abort does not reach the addon", () => {
    const aborted: number[] = [];

    const host = fetchHostWith((request, _scope, onReply) => {
      onReply({ kind: "reset", seq: request.seq, message: "When reading the response body: disk gone" });

      return {
        abort: () => {
          aborted.push(request.seq);
        },
      };
    });

    host.send({ kind: "request", request: REQUEST });
    host.send({ kind: "abort", seq: REQUEST.seq });

    expect([aborted, host.replies]).toStrictEqual([
      [],
      [{ kind: "reset", seq: 3, message: "When reading the response body: disk gone" }],
    ]);
  });
});

describe("backend_host with the renderer fetch over a message channel", () => {
  test("a renderer request reaches the addon and receives its response and its typed errors", async () => {
    const channel = new MessageChannel();
    const encoder = new TextEncoder();

    const fetch: AddonBackend["fetch"] = (request, _scope, onReply) => {
      const found = request.url === "/api/version";
      const body = found ? '{"version":"2.0.0"}' : '{"code":"not_found","message":"no run `r9`","causes":[]}';

      onReply({
        kind: "head",
        seq: request.seq,
        status: found ? 200 : 404,
        headers: [{ name: "content-type", value: "application/json" }],
      });
      onReply({ kind: "chunk", seq: request.seq, data: encoder.encode(body) });
      onReply({ kind: "end", seq: request.seq });

      return { abort: () => undefined };
    };

    serveFetch(nodeEndpoint(channel.port1), addonWith({ fetch }), "renderer");

    const portFetch = createPortFetch({
      onPort: (listener) => {
        listener(nodeFetchPort(channel.port2));
      },
      onStopped: () => undefined,
    });

    try {
      const version = await portFetch("http://treetime.desktop/api/version");
      const missing = await portFetch("http://treetime.desktop/api/runs/r9");

      expect([version.status, await version.json(), missing.status, await missing.json()]).toStrictEqual([
        200,
        { version: "2.0.0" },
        404,
        { code: "not_found", message: "no run `r9`", causes: [] },
      ]);
    } finally {
      channel.port1.close();
      channel.port2.close();
    }
  });
});

const REQUEST: PortRequest = { seq: 3, method: "GET", url: "/api/version", headers: [] };

function fetchHostWith(fetch: AddonBackend["fetch"], rejectRequest: AddonBackend["rejectRequest"] = rejectNothing) {
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

  serveFetch(endpoint, addonWith({ fetch, rejectRequest }), "renderer");

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

function rejectNothing(): PortReply[] {
  throw new Error("rejectRequest is not part of this test");
}

function addonWith(overrides: Partial<AddonBackend>): AddonBackend {
  return {
    rejectRequest: () => {
      throw new Error("rejectRequest is not part of this test");
    },
    fetch: () => {
      throw new Error("fetch is not part of this test");
    },
    ...overrides,
  };
}

type NodePort = InstanceType<typeof MessageChannel>["port1"];

function nodeEndpoint(port: NodePort): FetchEndpoint {
  return {
    post: (reply) => {
      port.postMessage(reply);
    },
    listen: (listener) => {
      port.on("message", (message: PortMessage) => {
        listener(message);
      });
      port.start();
    },
    onClose: (listener) => {
      port.on("close", listener);
    },
  };
}

function nodeFetchPort(port: NodePort): FetchPort {
  return {
    postMessage: (message) => {
      port.postMessage(message);
    },
    addEventListener: (_type, listener) => {
      port.on("message", (data: PortReply) => {
        listener({ data });
      });
    },
    start: () => {
      port.start();
    },
  };
}
