import { describe, expect, test } from "vitest";

import {
  portEndpoint,
  zBackendReply,
  zBackendRequest,
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
