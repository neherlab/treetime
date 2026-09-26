import {
  CancelledError,
  createBridge,
  bridgeErrorFromText,
  type BridgeTransport,
  type OperationRequestInput,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import { type ApiClient, runsUploadInput } from "@neherlab/app-contracts/client";
import { EventSourceParserStream } from "eventsource-parser/stream";

import { downloadBlob } from "./save-web";

interface WebBridgeDeps {
  client: ApiClient;
  fetchFn?: typeof fetch;
  apiBase?: string;
  saveBlob?: (blob: Blob, name: string) => void;
}

type Method = "GET" | "POST" | "PATCH" | "PUT" | "DELETE";

export function createWebBridge(deps: WebBridgeDeps): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const apiBase = deps.apiBase ?? "/api";
  const saveBlob = deps.saveBlob ?? downloadBlob;

  const debug = readDebugFlag();

  async function send(
    method: Method,
    path: string,
    body?: BodyInit,
    headers: HeadersInit = {},
    context = `${method} ${path}`,
  ): Promise<Response> {
    if (debug) console.debug("[TreeTime]", context);

    const init: RequestInit = { method, headers };

    if (body !== undefined) {
      init.body = body;
    }

    const response = await fetchFn(`${apiBase}/${path}`, init);

    if (!response.ok) {
      throw await responseError(response, `${context}: ${response.status}`);
    }

    return response;
  }

  async function call(request: OperationRequestInput): Promise<unknown> {
    const response = await send(
      "POST",
      "operations",
      JSON.stringify(request),
      { "Content-Type": "application/json" },
      request.operation,
    );

    const data: unknown = await response.json();

    if (debug) console.debug("[TreeTime]", request.operation, JSON.stringify(data));

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
    call,
    runEvents,
    saveRunFile: (id, path, name) => save(filePath(id, path), name),
    saveRunArchive: (id, name) => save(`${runPath(id)}/archive`, name),
    uploadInput: async (id, name, data) => {
      const { data: uploaded } = await runsUploadInput({
        client: deps.client,
        path: { id, name },
        body: data,
        throwOnError: true,
      });

      return uploaded;
    },
  };

  return createBridge(transport);
}

function filePath(id: string, path: string): string {
  return `${runPath(id)}/file?path=${encodeURIComponent(path)}`;
}

function runPath(id: string): string {
  return `runs/${encodeURIComponent(id)}`;
}

async function responseError(response: Response, context: string): Promise<Error> {
  const text = await response.text();

  return bridgeErrorFromText(text === "" ? response.statusText : text, context);
}

function readDebugFlag(): boolean {
  const env = import.meta.env;

  return env.VITE_TREETIME_DEBUG_FETCH === "true" || (env.DEV && env.VITE_TREETIME_DEBUG_FETCH !== "false");
}
