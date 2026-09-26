import { Backend } from "@neherlab/app-napi";
import type { MessagePortMain } from "electron";

import { saveRunFiles, serveBackend } from "./backend-host";
import { zBackendRequest, zControlRequest, type HostEndpoint, type SaveRequest } from "./backend-protocol";
import { DIAGNOSTIC_DIR_ENV, initDiagnostics } from "./diagnostics";

const diagnosticDir = process.env[DIAGNOSTIC_DIR_ENV];

if (diagnosticDir !== undefined) {
  initDiagnostics("treetime-backend", diagnosticDir);
}

const runsDir = process.argv.at(-1);

if (runsDir === undefined) {
  throw new Error("the back end needs the runs directory as its last argument");
}

const backend = new Backend(runsDir);

process.parentPort.on("message", (message) => {
  const control = zControlRequest.safeParse(message.data);

  if (!control.success) {
    console.warn("[TreeTime back end] ignored a malformed control message", control.error.message);

    return;
  }

  if (control.data.kind === "port") {
    const [port] = message.ports;

    if (port !== undefined) {
      serveBackend(portEndpoint(port), backend);
    }

    return;
  }

  void save(control.data);
});

async function save(request: SaveRequest): Promise<void> {
  const reply = await saveRunFiles(backend, request);
  // oxlint-disable-next-line unicorn/require-post-message-target-origin -- the parent port of a utility process takes no target origin
  process.parentPort.postMessage(reply);
}

function portEndpoint(port: MessagePortMain): HostEndpoint {
  return {
    post: (message) => {
      port.postMessage(message);
    },
    listen: (listener) => {
      port.on("message", (event) => {
        const request = zBackendRequest.safeParse(event.data);

        if (request.success) {
          listener(request.data);
        } else {
          console.warn("[TreeTime back end] ignored a malformed request", request.error.message);
        }
      });
      port.start();
    },
    onClose: (listener) => {
      port.on("close", listener);
    },
  };
}
