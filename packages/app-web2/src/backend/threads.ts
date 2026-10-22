import { errorMessage } from "@neherlab/app-contracts";

export const DEFAULT_JOBS = 8;

const POSITIVE_INTEGER = /^[1-9]\d*$/u;

const SHARED_MEMORY_PROBE: WebAssembly.MemoryDescriptor = { initial: 1, maximum: 65_536, shared: true };

export function parseJobs(search: string): JobsParameter {
  const value = new URLSearchParams(search).get("jobs");

  if (value === null) {
    return { jobs: DEFAULT_JOBS, warning: undefined };
  }

  const jobs = Number(value);

  if (POSITIVE_INTEGER.test(value) && Number.isSafeInteger(jobs)) {
    return { jobs, warning: undefined };
  }

  return {
    jobs: DEFAULT_JOBS,
    warning: `the URL parameter jobs=${value} is not a positive integer; runs use ${String(DEFAULT_JOBS)} threads`,
  };
}

export function chooseRunThreads(jobs: number, scope: ThreadScope, threadedWasmUrl: string): RunThreads {
  if (jobs === 1) {
    return { jobs, threadedWasmUrl: undefined, warning: undefined };
  }

  const missing = missingThreadSupport(scope);

  if (missing !== undefined) {
    return { jobs, threadedWasmUrl: undefined, warning: `runs use one thread, because ${missing}` };
  }

  return { jobs, threadedWasmUrl, warning: undefined };
}

export interface JobsParameter {
  jobs: number;
  warning: string | undefined;
}

export interface RunThreads {
  jobs: number;
  threadedWasmUrl: string | undefined;
  warning: string | undefined;
}

export interface ThreadScope {
  crossOriginIsolated: boolean;
  SharedArrayBuffer?: SharedArrayBufferConstructor | undefined;
  WebAssembly: { Memory: new (descriptor: WebAssembly.MemoryDescriptor) => { buffer: ArrayBufferLike } };
}

function missingThreadSupport(scope: ThreadScope): string | undefined {
  if (!scope.crossOriginIsolated) {
    return "the page is not cross-origin isolated";
  }

  if (scope.SharedArrayBuffer === undefined) {
    return "the browser has no SharedArrayBuffer";
  }

  try {
    const { buffer } = new scope.WebAssembly.Memory(SHARED_MEMORY_PROBE);

    return buffer instanceof scope.SharedArrayBuffer
      ? undefined
      : "the browser gives WebAssembly memory no shared buffer";
  } catch (error) {
    return `the browser cannot create shared WebAssembly memory: ${errorMessage(error)}`;
  }
}
