import {
  CancelledError,
  createBridge,
  parseJobEvent,
  type AppCommand,
  type BridgeTransport,
  type LogEvent,
  type TransportCommandOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import { EventSourceParserStream } from "eventsource-parser/stream";

export interface WebBridgeDeps {
  fetchFn?: typeof fetch;
  apiBase?: string;
}

export function createWebBridge(deps: WebBridgeDeps = {}): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const apiBase = deps.apiBase ?? "/api";

  const debug = readDebugFlag();

  async function getJson(path: string): Promise<unknown> {
    if (debug) console.debug("[TreeTime] GET", path);
    const response = await fetchFn(`${apiBase}/${path}`);

    if (!response.ok) {
      throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
    }

    const data: unknown = await response.json();

    if (debug) console.debug("[TreeTime] GET", path, JSON.stringify(data));

    return data;
  }

  async function postJson(path: string, body: unknown): Promise<unknown> {
    if (debug) console.debug("[TreeTime] POST", path);

    const response = await fetchFn(`${apiBase}/${path}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });

    if (!response.ok) {
      throw new Error(`POST ${path}: ${response.status} ${response.statusText}`);
    }

    const data: unknown = await response.json();

    return data;
  }

  async function postSse(command: AppCommand, config: unknown, options: TransportCommandOptions): Promise<unknown> {
    const { signal } = options;
    let jobId: string | undefined;

    const requestCancel = () => {
      if (jobId !== undefined) {
        void postJson(`jobs/${jobId}/cancel`, {}).catch((error: unknown) => {
          console.warn("[TreeTime] cancellation request failed", error);
        });
      }
    };

    signal?.addEventListener("abort", requestCancel);

    try {
      const response = await fetchFn(`${apiBase}/${command}`, {
        method: "POST",
        headers: { "Content-Type": "application/json", Accept: "text/event-stream" },
        body: JSON.stringify(config),
      });

      if (!response.ok || response.body === null) {
        throw new Error(`POST ${command}: ${response.status} ${response.statusText}`);
      }

      const messages = response.body.pipeThrough(new TextDecoderStream()).pipeThrough(new EventSourceParserStream());

      for await (const message of messages) {
        if (debug) console.debug("[TreeTime]", JSON.stringify(message));

        const data: unknown = JSON.parse(message.data);
        const event = parseJobEvent({ type: message.event, data });

        if (event.type === "terminal") {
          return event.data;
        }

        if (event.type === "started") {
          jobId = event.data.job_id;

          if (signal?.aborted === true) {
            requestCancel();
          }
        }

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
    } finally {
      signal?.removeEventListener("abort", requestCancel);
    }

    throw new Error(`${command}: the event stream ended without a terminal event`);
  }

  const transport: BridgeTransport = {
    query: (endpoint) => getJson(endpoint),
    request: (endpoint, body) => postJson(endpoint, body),
    command: (command, config, options) => postSse(command, config, options),
  };

  return createBridge(transport);
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
