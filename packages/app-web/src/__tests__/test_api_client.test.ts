import { afterEach, describe, expect, test, vi } from "vitest";

import { createWebApiClient } from "../api-client";

describe("api_client web", () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  test("the web client sends requests to the origin of the page", () => {
    vi.stubGlobal("location", { origin: "https://treetime.example" });

    expect(createWebApiClient().getConfig().baseUrl).toBe("https://treetime.example");
  });
});
