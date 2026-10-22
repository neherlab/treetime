import createBackendModule from "@neherlab/app-wasm";

import { JOB_WORKER_NAME } from "../job/protocol";
import { type BackendScope, serveBackend } from "./backend-worker";
import { compileWasm } from "./compile";

const scope: BackendScope = self;

serveBackend(
  scope,
  createBackendModule,
  () => new Worker(new URL("../job/worker.ts", import.meta.url), { type: "module", name: JOB_WORKER_NAME }),
  compileWasm,
);
