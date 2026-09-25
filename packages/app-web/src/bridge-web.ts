import {
  CancelledError,
  createBridge,
  bridgeErrorFromText,
  type BridgeTransport,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import { EventSourceParserStream } from "eventsource-parser/stream";

interface WebBridgeDeps {
  fetchFn?: typeof fetch;
  apiBase?: string;
  saveBlob?: (blob: Blob, name: string) => void;
}

type Method = "GET" | "POST" | "PATCH" | "PUT" | "DELETE";

export function createWebBridge(deps: WebBridgeDeps = {}): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const apiBase = deps.apiBase ?? "/api";
  const saveBlob = deps.saveBlob ?? downloadBlob;

  const debug = readDebugFlag();

  async function send(method: Method, path: string, body?: BodyInit, headers: HeadersInit = {}): Promise<Response> {
    if (debug) console.debug("[TreeTime]", method, path);

    const init: RequestInit = { method, headers };

    if (body !== undefined) {
      init.body = body;
    }

    const response = await fetchFn(`${apiBase}/${path}`, init);

    if (!response.ok) {
      throw await responseError(response, `${method} ${path}: ${response.status}`);
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

  async function save(path: string, name: string): Promise<boolean> {
    const response = await send("GET", path);
    saveBlob(await response.blob(), name);

    return true;
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
        throw await responseError(response, `GET ${path}: ${response.status}`);
      }

      const messages = response.body.pipeThrough(new TextDecoderStream()).pipeThrough(new EventSourceParserStream());

      for await (const message of messages) {
        if (debug) console.debug("[TreeTime]", JSON.stringify(message));

        const event: unknown = JSON.parse(message.data);

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
    saveRunFile: (id, path, name) => save(filePath(id, path), name),
    saveRunArchive: (id, name) => save(`${runPath(id)}/archive`, name),
    uploadInput: async (id, name, data) => {
      const response = await send("PUT", `${runPath(id)}/inputs/${encodeURIComponent(name)}`, data, {
        "Content-Type": "application/octet-stream",
      });

      const uploaded: unknown = await response.json();

      return uploaded;
    },
    runResults: (id) => json("GET", `${runPath(id)}/results`),
    runAuspice: (id) => json("GET", `${runPath(id)}/auspice`),
    compareRuns: (id, other) => json("GET", `${runPath(id)}/compare/${encodeURIComponent(other)}`),
    cladeInRuns: (request) => json("POST", "clade-in-runs", request),
  };

  return createBridge(transport);
}

function filePath(id: string, path: string): string {
  return `${runPath(id)}/file?path=${encodeURIComponent(path)}`;
}

function runPath(id: string): string {
  return `runs/${encodeURIComponent(id)}`;
}

function downloadBlob(blob: Blob, name: string): void {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");

  anchor.href = url;
  anchor.download = name;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
}

async function responseError(response: Response, context: string): Promise<Error> {
  const text = await response.text();

  return bridgeErrorFromText(text === "" ? response.statusText : text, context);
}

function readDebugFlag(): boolean {
  const env = import.meta.env;

  return env.VITE_TREETIME_DEBUG_FETCH === "true" || (env.DEV && env.VITE_TREETIME_DEBUG_FETCH !== "false");
}
