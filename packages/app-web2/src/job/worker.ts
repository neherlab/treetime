import createModule from "@neherlab/app-wasm";

import type { ModuleVariant } from "./protocol";
import { type JobModuleFactory, type JobScope, serveJob } from "./serve";

const scope: JobScope = self;

serveJob(scope, load);

async function load(variant: ModuleVariant, jobs: number): Promise<JobModuleFactory> {
  if (variant === "single") {
    return createModule;
  }

  const { default: createThreadedModule } = await import("@neherlab/app-wasm/threads");

  return (options) => {
    Object.assign(globalThis, { treetimeJobs: jobs });

    return createThreadedModule(options);
  };
}
