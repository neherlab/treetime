import { createWebApiClient } from "../api-client";

const RUNTIME_LOCATION = Object.getOwnPropertyDescriptor(globalThis, "location");

describe("api_client web", () => {
  afterEach(() => {
    if (RUNTIME_LOCATION === undefined) {
      Reflect.deleteProperty(globalThis, "location");
    } else {
      Object.defineProperty(globalThis, "location", RUNTIME_LOCATION);
    }
  });

  test("the web client sends requests to the origin of the page", () => {
    Object.defineProperty(globalThis, "location", {
      value: { origin: "https://treetime.example" },
      configurable: true,
    });

    expect(createWebApiClient().getConfig().baseUrl).toBe("https://treetime.example");
  });
});
