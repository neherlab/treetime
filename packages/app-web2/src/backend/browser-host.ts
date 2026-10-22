import wasmUrl from "@neherlab/app-wasm/app_wasm.wasm?url";
import threadedWasmUrl from "@neherlab/app-wasm/threads/app_wasm.wasm?url";

import { compileWasm } from "./compile";
import type { WorkerHost } from "./connection";
import { BACKEND_WORKER_NAME } from "./protocol";
import { chooseRunThreads } from "./threads";

export function browserWorkerHost(jobs: number): WorkerHost {
  const { threadedWasmUrl: threaded, warning } = chooseRunThreads(
    jobs,
    globalThis,
    new URL(threadedWasmUrl, globalThis.location.href).href,
  );

  if (warning !== undefined) {
    console.warn(`[TreeTime] ${warning}`);
  }

  return {
    start: () => new Worker(new URL("./worker.ts", import.meta.url), { type: "module", name: BACKEND_WORKER_NAME }),
    compile: () => compileWasm(wasmUrl),
    channel: () => new MessageChannel(),
    runThreads: { jobs, threadedWasmUrl: threaded },
  };
}
