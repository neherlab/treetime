import { createApiClient, type ApiClient } from "@neherlab/app-contracts/client";

const BASE_URL = "http://treetime.test";

export const RECORD = {
  id: "r1",
  title: "clock",
  command: "clock",
  config: { tree: "t.nwk" },
  status: "running",
  pinned: false,
  created_at: "2026-09-25T10:00:00Z",
  started_at: "2026-09-25T10:00:01Z",
  treetime_version: "1.0.0",
  inputs: [],
  changed_settings: [],
  headline: {},
  output_files: [],
} as const;

export interface Sent {
  key: string;
  body: string;
}

export type Route = (request: Request, attempt: number) => Response | Promise<Response>;

export class FakeServer {
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

    this.sent.push({ key, body: await request.text() });

    const route = this.routes[key];

    if (route === undefined) {
      return json({ code: "not_found", message: `no route for ${key}`, causes: [] }, 404);
    }

    return route(request, attempt);
  };

  client(): ApiClient {
    return createApiClient({ baseUrl: BASE_URL, fetch: this.fetch });
  }

  keys(): string[] {
    return this.sent.map((sent) => sent.key);
  }
}

export function json(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

function sseEvent(event: { seq: number }): string {
  return `id: ${event.seq}\ndata: ${JSON.stringify(event)}\n\n`;
}

export function eventStream(events: ReadonlyArray<{ seq: number }>, failure?: Error): Response {
  const encoder = new TextEncoder();
  const chunks = events.map(sseEvent);
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

export async function noDelay(): Promise<void> {
  await Promise.resolve();
}
