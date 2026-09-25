import { Backend } from "@neherlab/app-napi";
import type { MessagePortMain } from "electron";

import { serveBackend } from "./backend-host";
import { zBackendRequest, type HostEndpoint } from "./backend-protocol";
import { initDiagnostics } from "./diagnostics";

initDiagnostics("treetime-backend");

const runsDir = process.argv.at(-1);

if (runsDir === undefined) {
  throw new Error("the back end needs the runs directory as its last argument");
}

const backend = new Backend(runsDir);

process.parentPort.on("message", (message) => {
  const [port] = message.ports;

  if (port !== undefined) {
    serveBackend(portEndpoint(port), backend);
  }
});

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
