import { type EmscriptenOptions, instantiateBackend } from "../instantiate";

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
    const exported = await instantiateBackend(MODULE_EXPORTING_SEVEN, receivedExport);

    expect(exported).toBe(7);
  });

  test("a runtime that cannot take the instance rejects the instantiation", async () => {
    const instantiation = instantiateBackend(MODULE_EXPORTING_SEVEN, rejectingRuntime);

    await expect(instantiation).rejects.toThrow("the module has no memory");
  });

  test("a module whose imports the runtime does not provide rejects the instantiation", async () => {
    const instantiation = instantiateBackend(MODULE_IMPORTING_A_FUNCTION, receivedExport);

    await expect(instantiation).rejects.toBeInstanceOf(TypeError);
  });
});

function receivedExport({ instantiateWasm }: EmscriptenOptions): Promise<unknown> {
  return new Promise((resolve) => {
    instantiateWasm({}, (instance) => {
      const global = instance.exports["g"];

      resolve(global instanceof WebAssembly.Global ? global.value : undefined);
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
