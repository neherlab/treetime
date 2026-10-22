import { type EmscriptenOptions, ErrorLines, instantiateBackend } from "../instantiate";

const MODULE_EXPORTING_SEVEN = new WebAssembly.Module(
  new Uint8Array([
    0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00, 0x06, 0x06, 0x01, 0x7f, 0x00, 0x41, 0x07, 0x0b, 0x07, 0x05, 0x01,
    0x01, 0x67, 0x03, 0x00,
  ]),
);

const MODULE_IMPORTING_A_FUNCTION = new WebAssembly.Module(
  new Uint8Array([
    0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00, 0x01, 0x04, 0x01, 0x60, 0x00, 0x00, 0x02, 0x09, 0x01, 0x03, 0x65,
    0x6e, 0x76, 0x01, 0x66, 0x00, 0x00,
  ]),
);

describe("backend instantiation", () => {
  test("the runtime receives an instance of the compiled module", async () => {
    const exported = await instantiateBackend(MODULE_EXPORTING_SEVEN, receivedExport, errorLines());

    expect(exported).toBe(7);
  });

  test("the runtime receives the compiled module with the instance, which it passes on to its pthread workers", async () => {
    const received = await instantiateBackend(MODULE_EXPORTING_SEVEN, receivedModule, errorLines());

    expect(received).toBe(MODULE_EXPORTING_SEVEN);
  });

  test("a runtime that cannot take the instance rejects the instantiation", async () => {
    const instantiation = instantiateBackend(MODULE_EXPORTING_SEVEN, rejectingRuntime, errorLines());

    await expect(instantiation).rejects.toThrow("the module has no memory");
  });

  test("a module whose imports the runtime does not provide rejects the instantiation", async () => {
    const instantiation = instantiateBackend(MODULE_IMPORTING_A_FUNCTION, receivedExport, errorLines());

    await expect(instantiation).rejects.toBeInstanceOf(TypeError);
  });
});

describe("backend error lines", () => {
  test("the runtime's stderr lines, the abort reason, and the error message form the error lines in this order", async () => {
    const errors = errorLines();

    errors.error("RuntimeError: unreachable");
    await instantiateBackend(
      MODULE_EXPORTING_SEVEN,
      (options) => {
        options.printErr("memory allocation of 4294967296 bytes failed");
        options.onAbort("stack overflow");

        return receivedExport(options);
      },
      errors,
    );

    expect(errors.lines()).toStrictEqual([
      "memory allocation of 4294967296 bytes failed",
      "stack overflow",
      "RuntimeError: unreachable",
    ]);
  });

  test("the error lines keep the last 20 stderr lines", () => {
    const errors = errorLines();

    Array.from({ length: 25 }, (_, index) => `line ${String(index + 1)}`).forEach((line) => {
      errors.printErr(line);
    });

    expect([errors.lines().length, errors.lines()[0], errors.lines().at(-1)]).toStrictEqual([20, "line 6", "line 25"]);
  });

  test("an abort without a reason adds no line and still reports the abort", () => {
    let aborts = 0;

    const errors = new ErrorLines(() => {
      aborts += 1;
    });

    errors.abort(undefined);

    expect([errors.lines(), aborts]).toStrictEqual([[], 1]);
  });
});

function errorLines(): ErrorLines {
  return new ErrorLines(() => undefined);
}

function receivedExport({ instantiateWasm }: EmscriptenOptions): Promise<unknown> {
  return new Promise((resolve) => {
    instantiateWasm({}, (instance) => {
      const global = instance.exports["g"];

      resolve(global instanceof WebAssembly.Global ? global.value : undefined);
    });
  });
}

function receivedModule({ instantiateWasm }: EmscriptenOptions): Promise<WebAssembly.Module> {
  return new Promise((resolve) => {
    instantiateWasm({}, (_instance, module) => {
      resolve(module);
    });
  });
}

function rejectingRuntime({ instantiateWasm }: EmscriptenOptions): Promise<never> {
  return new Promise(() => {
    instantiateWasm({}, () => {
      throw new Error("the module has no memory");
    });
  });
}
