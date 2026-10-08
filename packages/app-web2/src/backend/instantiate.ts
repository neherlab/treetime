export function instantiateBackend<T>(
  module: WebAssembly.Module,
  create: (options: EmscriptenOptions) => Promise<T>,
): Promise<T> {
  const failed = Promise.withResolvers<never>();

  const options: EmscriptenOptions = {
    instantiateWasm: (imports, receive) => {
      WebAssembly.instantiate(module, imports).then(receive).catch(failed.reject);
    },
  };

  return Promise.race([create(options), failed.promise]);
}

export interface EmscriptenOptions {
  instantiateWasm: (imports: WebAssembly.Imports, receive: (instance: WebAssembly.Instance) => void) => void;
}
