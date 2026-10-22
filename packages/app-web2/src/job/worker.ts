import createModule from "@neherlab/app-wasm";

import type { ModuleVariant } from "./protocol";
import { type JobModuleFactory, type JobScope, serveJob } from "./serve";

const scope: JobScope = self;

serveJob(scope, load);

function load(variant: ModuleVariant): Promise<JobModuleFactory> {
  if (variant === "threaded") {
    return Promise.reject(new Error("this build has no threaded WebAssembly module"));
  }

  return Promise.resolve(createModule);
}
