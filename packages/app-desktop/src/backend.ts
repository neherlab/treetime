import { Backend, type PortMessage, type PortReply } from "@neherlab/app-napi";
import type { MessagePortMain } from "electron";

import { serveFetch } from "./backend-host";
import type { ControlReply, ControlRequest, FetchEndpoint } from "./backend-protocol";
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

  const [port] = message.ports;

  if (port !== undefined) {
    serveFetch(fetchEndpoint(port), backend, control.scope);
  }
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
