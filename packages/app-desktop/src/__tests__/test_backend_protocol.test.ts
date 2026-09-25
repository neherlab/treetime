import { describe, expect, test } from "vitest";

import {
  portEndpoint,
  zBackendReply,
  zBackendRequest,
  zControlReply,
  zControlRequest,
  type BackendReply,
  type BackendRequest,
  type PortLike,
} from "../backend-protocol";

describe("backend_protocol requests", () => {
  test.each([
    [{ kind: "call", seq: 1, request: "{}" }],
    [{ kind: "subscribe", seq: 2, id: "r1", from: 3 }],
    [{ kind: "unsubscribe", seq: 2 }],
  ])("%j is a request", (request) => {
    expect(zBackendRequest.parse({ ...request, extra: true })).toStrictEqual(request);
  });

  test.each([
    [null],
    ["call"],
    [{ kind: "call", seq: "1", request: "{}" }],
    [{ kind: "call", seq: -1, request: "{}" }],
    [{ kind: "call", seq: 1, request: {} }],
    [{ kind: "subscribe", seq: 1, id: "r1", from: "0" }],
    [{ kind: "format-disk", seq: 1 }],
  ])("%j is not a request", (data) => {
    expect(zBackendRequest.safeParse(data).success).toBe(false);
  });
});

describe("backend_protocol replies", () => {
  test.each([
    [{ kind: "result", seq: 1, json: "{}" }],
    [{ kind: "event", seq: 1, json: "{}" }],
    [{ kind: "error", seq: 1, error: "failed" }],
  ])("%j is a reply", (reply) => {
    expect(zBackendReply.parse(reply)).toStrictEqual(reply);
  });

  test.each([
    [undefined],
    [{ kind: "result", seq: 1 }],
    [{ kind: "error", seq: 1, error: 1 }],
    [{ kind: "call", seq: 1, request: "{}" }],
  ])("%j is not a reply", (data) => {
    expect(zBackendReply.safeParse(data).success).toBe(false);
  });
});

describe("backend_protocol port endpoints", () => {
  test("a port endpoint delivers valid messages and drops malformed ones", () => {
    const listeners: Array<(event: { data: unknown }) => void> = [];

    const port: PortLike = {
      postMessage: () => undefined,
      addEventListener: (_type, listener) => {
        listeners.push(listener);
      },
      start: () => undefined,
    };

    const received: BackendReply[] = [];
    const endpoint = portEndpoint<BackendReply, BackendRequest>(port, zBackendReply);

    endpoint.listen((reply) => {
      received.push(reply);
    });
    listeners.forEach((listener) => {
      listener({ data: { kind: "result", seq: "x" } });
      listener({ data: { kind: "result", seq: 3, json: "{}" } });
    });

    expect(received).toStrictEqual([{ kind: "result", seq: 3, json: "{}" }]);
  });
});

describe("backend_protocol control messages", () => {
  test.each([
    [{ kind: "port" }],
    [{ kind: "save-file", seq: 1, id: "r1", path: "a.nwk", destination: "/home/user/a.nwk" }],
    [{ kind: "save-archive", seq: 2, id: "r1", destination: "/home/user/r1.zip" }],
  ])("%j is a control request", (request) => {
    expect(zControlRequest.parse(request)).toStrictEqual(request);
  });

  test.each([[{ kind: "saved", seq: 1 }], [{ kind: "error", seq: 1, error: "disk full" }]])(
    "%j is a control reply",
    (reply) => {
      expect(zControlReply.parse(reply)).toStrictEqual(reply);
    },
  );

  test("a save without its destination is not a control request", () => {
    expect(zControlRequest.safeParse({ kind: "save-archive", seq: 2, id: "r1" }).success).toBe(false);
  });
});
