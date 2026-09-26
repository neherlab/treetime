import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import {
  ApiError,
  createApiClient,
  runsCreate,
  runsEvents,
  runsGet,
  runsUploadInput,
  SSE_MAX_RETRY_ATTEMPTS,
  version,
} from "../client";
import type { RunEvent } from "../generated/types.gen";
import { zRunsUploadInputBody } from "../generated/zod.gen";

const BASE_URL = "http://treetime.test";

const NOT_FOUND = { code: "not_found", message: "When reading run `r9`", causes: ["no run with id `r9`"] };

const RECORD = {
  id: "r1",
  title: "clock",
  command: "clock",
  config: { tree: "t.nwk" },
  status: "created",
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

interface Sent {
  key: string;
  method: string;
  url: string;
  headers: Headers;
  body: Blob;
}

type Route = (request: Request, attempt: number) => Response;

class FakeServer {
  readonly sent: Sent[] = [];
  private readonly routes: Record<string, Route>;

  constructor(routes: Record<string, Route>) {
    this.routes = routes;
  }

  readonly fetch = async (input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
    const request = new Request(input, init);
    const url = new URL(request.url);
    const key = `${request.method} ${url.pathname}${url.search}`;
    const attempt = this.sent.filter((sent) => sent.key === key).length;

    this.sent.push({
      key,
      method: request.method,
      url: request.url,
      headers: request.headers,
      body: await request.blob(),
    });

    const route = this.routes[key];

    if (route === undefined) {
      return json({ code: "not_found", message: `no route for ${key}`, causes: [] }, 404);
    }

    return route(request, attempt);
  };
}

function json(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

function logEvent(seq: number, message: string): string {
  const data = { seq, time: "t", type: "log", data: { level: "info", message } };

  return `id: ${seq}\nevent: log\ndata: ${JSON.stringify(data)}\n\n`;
}

function eventStream(chunks: string[], failure?: Error): Response {
  const encoder = new TextEncoder();
  let next = 0;

  const body = new ReadableStream<Uint8Array>({
    pull(controller) {
      const chunk = chunks[next];
      next += 1;

      if (chunk !== undefined) {
        controller.enqueue(encoder.encode(chunk));
      } else if (failure === undefined) {
        controller.close();
      } else {
        controller.error(failure);
      }
    },
  });

  return new Response(body, { headers: { "Content-Type": "text/event-stream" } });
}

async function collect<T>(stream: AsyncGenerator<T>): Promise<T[]> {
  const items: T[] = [];

  for await (const item of stream) {
    items.push(item);
  }

  return items;
}

function messages(events: RunEvent[]): string[] {
  return events.flatMap((event) => (event.type === "log" ? [event.data.message] : []));
}

describe("api_client requests", () => {
  test("request URLs are the OpenAPI paths below the base URL", async () => {
    const server = new FakeServer({ "GET /api/runs/a%20b": () => json({ ...RECORD, id: "a b" }) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { data } = await runsGet({ client, path: { id: "a b" }, throwOnError: true });

    expect(data.id).toBe("a b");
    expect(server.sent.map((sent) => sent.url)).toStrictEqual([`${BASE_URL}/api/runs/a%20b`]);
  });

  test("a request body is sent as JSON exactly as given", async () => {
    const server = new FakeServer({ "POST /api/runs": () => json(RECORD, 201) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });
    const body = { command: "clock" as const, config: { tree: "t.nwk" } };

    const { data } = await runsCreate({ client, body, throwOnError: true });

    expect(data).toStrictEqual(RECORD);
    expect(server.sent[0]?.headers.get("Content-Type")).toBe("application/json");
    expect(JSON.parse(await (server.sent[0]?.body ?? new Blob()).text())).toStrictEqual(body);
  });

  test("an unknown field in a closed request object is rejected before sending", async () => {
    const server = new FakeServer({ "POST /api/runs": () => json(RECORD, 201) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });
    const body = { command: "clock" as const, config: {}, extra: 1 };

    await expect(runsCreate({ client, body, throwOnError: true })).rejects.toBeInstanceOf(ZodError);
    expect(server.sent).toStrictEqual([]);
  });

  test("a malformed response rejects with a ZodError", async () => {
    const server = new FakeServer({ "GET /api/version": () => json({ version: 1 }) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    await expect(version({ client, throwOnError: true })).rejects.toBeInstanceOf(ZodError);
  });

  test("an error answer rejects with an ApiError carrying the status and the ErrorResponse", async () => {
    const server = new FakeServer({ "GET /api/runs/r9": () => json(NOT_FOUND, 404) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const error = await runsGet({ client, path: { id: "r9" }, throwOnError: true }).catch(
      (failure: unknown) => failure,
    );

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({
      status: 404,
      response: NOT_FOUND,
      message: "When reading run `r9`: no run with id `r9`",
    });
  });

  test("an error answer without an ErrorResponse body becomes an internal error with the body text", async () => {
    const server = new FakeServer({ "GET /api/runs/r9": () => new Response("Bad Gateway", { status: 502 }) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const error = await runsGet({ client, path: { id: "r9" }, throwOnError: true }).catch(
      (failure: unknown) => failure,
    );

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({
      status: 502,
      response: { code: "internal_error", message: "Bad Gateway", causes: [] },
    });
  });

  test("without throwOnError, the error field holds the ApiError", async () => {
    const server = new FakeServer({ "GET /api/runs/r9": () => json(NOT_FOUND, 404) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const result = await runsGet({ client, path: { id: "r9" } });

    expect(result.data).toBeUndefined();
    expect(result.error).toBeInstanceOf(ApiError);
  });

  test("an upload sends the Blob bytes as application/octet-stream", async () => {
    const uploaded = { name: "a b.nwk", path: "/runs/r1/inputs/a b.nwk", size: 6, sha256: "x" };
    const server = new FakeServer({ "PUT /api/runs/r1/inputs/a%20b.nwk": () => json(uploaded) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { data } = await runsUploadInput({
      client,
      path: { id: "r1", name: "a b.nwk" },
      body: new Blob(["(A,B);"]),
      throwOnError: true,
    });

    expect(data).toStrictEqual(uploaded);
    expect(server.sent[0]?.headers.get("Content-Type")).toBe("application/octet-stream");
    expect(await (server.sent[0]?.body ?? new Blob()).text()).toBe("(A,B);");
  });

  test("an upload body validates as a Blob and nothing else", () => {
    expect(zRunsUploadInputBody.safeParse(new Blob(["(A,B);"])).success).toBe(true);
    expect(zRunsUploadInputBody.safeParse("(A,B);").success).toBe(false);
  });
});

describe("api_client event streams", () => {
  test("a dropped stream reconnects with Last-Event-ID and resumes after the last event", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events?from=0": (_request, attempt) =>
        attempt === 0
          ? eventStream([logEvent(0, "a"), logEvent(1, "b")], new TypeError("connection reset"))
          : eventStream([logEvent(2, "c")]),
    });

    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { stream } = await runsEvents({ client, path: { id: "r1" }, query: { from: 0 }, sseDefaultRetryDelay: 0 });

    expect(messages(await collect(stream))).toStrictEqual(["a", "b", "c"]);
    expect(server.sent.map((sent) => sent.headers.get("Last-Event-ID"))).toStrictEqual([null, "1"]);
  });

  test("an invalid event rejects the stream without reconnecting", async () => {
    const invalid = `id: 1\nevent: log\ndata: ${JSON.stringify({ seq: 1, time: "t", type: "log", data: { level: "loud" } })}\n\n`;
    const server = new FakeServer({ "GET /api/runs/r1/events": () => eventStream([logEvent(0, "a"), invalid]) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { stream } = await runsEvents({ client, path: { id: "r1" }, sseDefaultRetryDelay: 0 });
    const received: RunEvent[] = [];

    await expect(
      (async () => {
        for await (const event of stream) {
          received.push(event);
        }
      })(),
    ).rejects.toBeInstanceOf(ZodError);
    expect(messages(received)).toStrictEqual(["a"]);
    expect(server.sent).toHaveLength(1);
  });

  test("a client error answer rejects the stream with an ApiError without reconnecting", async () => {
    const server = new FakeServer({ "GET /api/runs/r9/events": () => json(NOT_FOUND, 404) });
    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { stream } = await runsEvents({ client, path: { id: "r9" }, sseDefaultRetryDelay: 0 });
    const error = await collect(stream).catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 404, response: NOT_FOUND });
    expect(server.sent).toHaveLength(1);
  });

  test("a server error answer is retried", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events": (_request, attempt) =>
        attempt === 0
          ? json({ code: "internal_error", message: "busy", causes: [] }, 503)
          : eventStream([logEvent(0, "a")]),
    });

    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { stream } = await runsEvents({ client, path: { id: "r1" }, sseDefaultRetryDelay: 0 });

    expect(messages(await collect(stream))).toStrictEqual(["a"]);
    expect(server.sent).toHaveLength(2);
  });

  test("connection attempts stop after the retry limit", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events": () => eventStream([], new TypeError("connection reset")),
    });

    const client = createApiClient({ baseUrl: BASE_URL, fetch: server.fetch });

    const { stream } = await runsEvents({ client, path: { id: "r1" }, sseDefaultRetryDelay: 0 });

    expect(await collect(stream)).toStrictEqual([]);
    expect(server.sent).toHaveLength(SSE_MAX_RETRY_ATTEMPTS);
  });
});
