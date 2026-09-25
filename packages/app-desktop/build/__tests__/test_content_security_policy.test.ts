import { describe, expect, test } from "vitest";

import { contentSecurityPolicy } from "../content-security-policy";

describe("content security policy", () => {
  test("the packaged application allows scripts of its own origin only", () => {
    expect(contentSecurityPolicy(false)).toBe(
      "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; " +
        "font-src 'self' data:; connect-src 'self'; worker-src 'self' blob:; object-src 'none'; base-uri 'none'; " +
        "form-action 'none'",
    );
  });

  test("the development server also allows the inline script of hot reloading", () => {
    expect(contentSecurityPolicy(true)).toContain("script-src 'self' 'unsafe-inline';");
  });
});
