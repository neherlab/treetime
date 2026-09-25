import {
  CancelledError,
  createBridge,
  parseRunEvent,
  zErrorResponse,
  type BridgeTransport,
  type LogEvent,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import { EventSourceParserStream } from "eventsource-parser/stream";

export interface WebBridgeDeps {
  fetchFn?: typeof fetch;
  apiBase?: string;
}

type Method = "GET" | "POST" | "PATCH" | "PUT" | "DELETE";

export function createWebBridge(deps: WebBridgeDeps = {}): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const apiBase = deps.apiBase ?? "/api";

  const debug = readDebugFlag();

  async function send(method: Method, path: string, body?: BodyInit, headers: HeadersInit = {}): Promise<Response> {
    if (debug) console.debug("[TreeTime]", method, path);

    const init: RequestInit = { method, headers };

    if (body !== undefined) {
      init.body = body;
    }

    const response = await fetchFn(`${apiBase}/${path}`, init);

    if (!response.ok) {
      throw new Error(`${method} ${path}: ${response.status} ${await errorMessage(response)}`);
    }

    return response;
  }

  async function json(method: Method, path: string, body?: unknown): Promise<unknown> {
    const response =
      body === undefined
        ? await send(method, path)
        : await send(method, path, JSON.stringify(body), { "Content-Type": "application/json" });

    const data: unknown = await response.json();

    if (debug) console.debug("[TreeTime]", method, path, JSON.stringify(data));

    return data;
  }

  async function bytes(path: string): Promise<Uint8Array> {
    const response = await send("GET", path);

    return new Uint8Array(await response.arrayBuffer());
  }

  async function runEvents(id: string, options: TransportEventOptions): Promise<void> {
    const init: RequestInit = { headers: { Accept: "text/event-stream" } };

    if (options.signal !== undefined) {
      init.signal = options.signal;
    }

    const path = `runs/${encodeURIComponent(id)}/events?from=${options.from}`;

    try {
      const response = await fetchFn(`${apiBase}/${path}`, init);

      if (!response.ok || response.body === null) {
        throw new Error(`GET ${path}: ${response.status} ${await errorMessage(response)}`);
      }

      const messages = response.body.pipeThrough(new TextDecoderStream()).pipeThrough(new EventSourceParserStream());

      for await (const message of messages) {
        if (debug) console.debug("[TreeTime]", JSON.stringify(message));

        const event = parseRunEvent(JSON.parse(message.data));

        if (event.type === "log") {
          logToConsole(event.data);
        }

        options.onEvent(event);
      }
    } catch (err: unknown) {
      if (err instanceof DOMException && err.name === "AbortError") {
        throw new CancelledError();
      }

      throw err;
    }
  }

  const transport: BridgeTransport = {
    version: () => json("GET", "version"),
    datasets: () => json("GET", "datasets"),
    checkConfig: (request) => json("POST", "check-config", request),
    runConfig: (request) => json("POST", "run-config", request),
    checkInputs: (request) => json("POST", "check-inputs", request),
    listRuns: () => json("GET", "runs"),
    createRun: (request) => json("POST", "runs", request),
    getRun: (id) => json("GET", runPath(id)),
    startRun: (id, request) => json("POST", `${runPath(id)}/start`, request),
    updateRun: (id, request) => json("PATCH", runPath(id), request),
    cancelRun: (id) => json("POST", `${runPath(id)}/cancel`),
    deleteRun: async (id) => {
      await send("DELETE", runPath(id));
    },
    restoreRun: (id) => json("POST", `${runPath(id)}/restore`),
    purgeRun: async (id) => {
      await send("POST", `${runPath(id)}/purge`);
    },
    runEvents,
    runFiles: (id) => json("GET", `${runPath(id)}/files`),
    readRunFile: (id, path) => bytes(`${runPath(id)}/file?path=${encodeURIComponent(path)}`),
    runArchive: (id) => bytes(`${runPath(id)}/archive`),
    uploadInput: async (id, name, data) => {
      const response = await send("PUT", `${runPath(id)}/inputs/${encodeURIComponent(name)}`, data, {
        "Content-Type": "application/octet-stream",
      });

      const uploaded: unknown = await response.json();

      return uploaded;
    },
  };

  return createBridge(transport);
}

function runPath(id: string): string {
  return `runs/${encodeURIComponent(id)}`;
}

async function errorMessage(response: Response): Promise<string> {
  const text = await response.text();
  const parsed = zErrorResponse.safeParse(parseJson(text));

  if (parsed.success) {
    return parsed.data.message;
  }

  return text === "" ? response.statusText : text;
}

function parseJson(text: string): unknown {
  try {
    const value: unknown = JSON.parse(text);

    return value;
  } catch {
    return undefined;
  }
}

function readDebugFlag(): boolean {
  const env = import.meta.env;

  return env.VITE_TREETIME_DEBUG_FETCH === "true" || (env.DEV && env.VITE_TREETIME_DEBUG_FETCH !== "false");
}

function logToConsole(log: LogEvent): void {
  switch (log.level) {
    case "error":
      console.error(`[TreeTime] ${log.message}`);
      break;
    case "warn":
      console.warn(`[TreeTime] ${log.message}`);
      break;
    case "info":
    case "debug":
    case "trace":
      console.log(`[TreeTime] [${log.level}] ${log.message}`);
      break;
  }
}
