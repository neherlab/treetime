import { Backend, type PortMessage, type PortReply } from "@neherlab/app-napi";
import type { MessagePortMain } from "electron";

import { saveRunFiles, serveFetch } from "./backend-host";
import type { ControlRequest, FetchEndpoint, SaveRequest } from "./backend-protocol";
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

process.parentPort.on("message", (message: { data: ControlRequest; ports: MessagePortMain[] }) => {
  const control = message.data;

  if (control.kind === "port") {
    const [port] = message.ports;

    if (port !== undefined) {
      serveFetch(fetchEndpoint(port), backend);
    }

    return;
  }

  void save(control);
});

async function save(request: SaveRequest): Promise<void> {
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
