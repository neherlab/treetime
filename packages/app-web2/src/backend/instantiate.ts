import { errorMessage } from "@neherlab/app-contracts";

const KEPT_ERROR_LINES = 20;

export function instantiateBackend<T>(
  module: WebAssembly.Module,
  create: (options: EmscriptenOptions) => Promise<T>,
  errors: ErrorLines,
): Promise<T> {
  const failed = Promise.withResolvers<never>();

  const options: EmscriptenOptions = {
    instantiateWasm: (imports, receive) => {
      WebAssembly.instantiate(module, imports).then(receive).catch(failed.reject);
    },
    printErr: (text) => {
      errors.printErr(text);
    },
    onAbort: (what) => {
      errors.abort(what);
    },
  };

  return Promise.race([create(options), failed.promise]);
}

export class ErrorLines {
  readonly #onAbort: () => void;
  readonly #stderr: string[] = [];
  #abort: string | undefined;
  #error: string | undefined;

  constructor(onAbort: () => void) {
    this.#onAbort = onAbort;
  }

  printErr(text: string): void {
    console.error(text);
    this.#stderr.push(text);
    this.#stderr.splice(0, Math.max(0, this.#stderr.length - KEPT_ERROR_LINES));
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- Emscripten passes the abort reason as any value, often a string
  abort(what: unknown): void {
    const reason = what === undefined ? "" : errorMessage(what);

    if (reason !== "") {
      this.#abort ??= reason;
    }

    this.#onAbort();
  }

  error(message: string): void {
    this.#error ??= message;
  }

  lines(): string[] {
    return [...this.#stderr, ...optional(this.#abort), ...optional(this.#error)];
  }
}

export interface EmscriptenOptions {
  instantiateWasm: (imports: WebAssembly.Imports, receive: (instance: WebAssembly.Instance) => void) => void;
  printErr: (text: string) => void;
  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- Emscripten passes the abort reason as any value
  onAbort: (what: unknown) => void;
}

function optional(text: string | undefined): string[] {
  return text === undefined ? [] : [text];
}
