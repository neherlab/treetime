import { describe, expect, test } from "vitest";

import { contentSecurityPolicy, isAppUrl, isExternalLink, isTrustedSender } from "../security";

const PACKAGED = "file:///opt/treetime/resources/app/dist/index.html";

const DEV_SERVER = "http://localhost:5173/";

describe("security app urls", () => {
  test.each([
    [PACKAGED, PACKAGED, true],
    [`${PACKAGED}#/runs/r1`, PACKAGED, true],
    ["file:///home/user/secrets.html", PACKAGED, false],
    ["http://localhost:5173/runs/r1", DEV_SERVER, true],
    ["http://localhost:5174/", DEV_SERVER, false],
    ["https://example.org/", DEV_SERVER, false],
    ["not a url", DEV_SERVER, false],
  ])("%s belongs to the application at %s: %s", (url, appUrl, expected) => {
    expect(isAppUrl(url, appUrl)).toBe(expected);
  });
});

describe("security senders", () => {
  test("the top-level frame of the application is trusted", () => {
    expect(isTrustedSender({ url: `${PACKAGED}#/`, parent: null }, PACKAGED)).toBe(true);
  });

  test("a subframe of the application is not trusted", () => {
    expect(isTrustedSender({ url: PACKAGED, parent: {} }, PACKAGED)).toBe(false);
  });

  test("a frame that navigated elsewhere is not trusted", () => {
    expect(isTrustedSender({ url: "https://example.org/", parent: null }, PACKAGED)).toBe(false);
  });

  test("a message without a frame is not trusted", () => {
    expect(isTrustedSender(null, PACKAGED)).toBe(false);
  });
});

describe("security external links", () => {
  test.each([
    ["https://doi.org/10.1093/ve/vex042", true],
    ["http://example.org", true],
    ["file:///etc/passwd", false],
    ["javascript:alert(1)", false],
  ])("%s opens in the system browser: %s", (url, expected) => {
    expect(isExternalLink(url)).toBe(expected);
  });
});

describe("security content policy", () => {
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
