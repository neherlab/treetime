import { Backend, type PortMessage, type PortReply } from "@neherlab/app-napi";
import type { MessagePortMain } from "electron";

import { saveRunFiles, serveFetch } from "./backend-host";
import type { ControlReply, ControlRequest, FetchEndpoint, SaveRequest } from "./backend-protocol";
import { DIAGNOSTIC_DIR_ENV, initDiagnostics } from "./diagnostics";
import { napiErrorResponse } from "./napi-error";

const diagnosticDir = process.env[DIAGNOSTIC_DIR_ENV];

if (diagnosticDir !== undefined) {
  initDiagnostics("treetime-backend", diagnosticDir);
}

const backend = startBackend();

process.parentPort.on("message", (message: { data: ControlRequest; ports: MessagePortMain[] }) => {
  const control = message.data;

  if (backend === undefined) {
    return;
  }

  if (control.kind === "port") {
    const [port] = message.ports;

    if (port !== undefined) {
      serveFetch(fetchEndpoint(port), backend);
    }

    return;
  }

  void save(backend, control);
});

function startBackend(): Backend | undefined {
  try {
    return new Backend();
  } catch (error: unknown) {
    // oxlint-disable-next-line unicorn/require-post-message-target-origin -- the parent port of a utility process takes no target origin
    process.parentPort.postMessage({ kind: "failed", error: napiErrorResponse(error) } satisfies ControlReply);

    return undefined;
  }
}

async function save(backend: Backend, request: SaveRequest): Promise<void> {
  const reply = await saveRunFiles(backend, request);
  // oxlint-disable-next-line unicorn/require-post-message-target-origin -- the parent port of a utility process takes no target origin
  process.parentPort.postMessage(reply);
}

function fetchEndpoint(port: MessagePortMain): FetchEndpoint {
  return {
    post: (reply: PortReply) => {
      port.postMessage(reply);
    },
    listen: (listener) => {
      port.on("message", (event: { data: PortMessage }) => {
        listener(event.data);
      });
      port.start();
    },
    onClose: (listener) => {
      port.on("close", listener);
    },
  };
}
