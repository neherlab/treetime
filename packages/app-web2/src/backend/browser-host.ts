import wasmUrl from "@neherlab/app-wasm/app_wasm.wasm?url";

import type { WorkerHost } from "./connection";
import { BACKEND_WORKER_NAME } from "./protocol";

const WASM_FETCH_TIMEOUT_MS = 120_000;

export const browserWorkerHost: WorkerHost = {
  start: () => new Worker(new URL("./worker.ts", import.meta.url), { type: "module", name: BACKEND_WORKER_NAME }),
  compile: async () => {
    const response = await fetch(wasmUrl, { signal: AbortSignal.timeout(WASM_FETCH_TIMEOUT_MS) });

    if (!response.ok) {
      throw new Error(`the WebAssembly module answered with HTTP status ${String(response.status)}`);
    }

    return WebAssembly.compileStreaming(response);
  },
  channel: () => new MessageChannel(),
};
