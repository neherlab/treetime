import { describe, expect, test } from "vitest";

import { APP_URL } from "../app-scheme";
import { confineNavigation, isTrustedSender } from "../security";

const PACKAGED = APP_URL;

const DEV_SERVER = "http://localhost:5173/";

describe("security navigation", () => {
  test.each([
    [PACKAGED, PACKAGED, true],
    ["app://treetime/runs/r1/results", PACKAGED, true],
    ["app://other/", PACKAGED, false],
    ["app:///index.html", PACKAGED, false],
    ["file:///opt/treetime/resources/app.asar/dist/index.html", PACKAGED, false],
    ["file:///home/user/secrets.html", PACKAGED, false],
    ["http://localhost:5173/runs/r1", DEV_SERVER, true],
    ["http://localhost:5174/", DEV_SERVER, false],
    ["https://example.org/", DEV_SERVER, false],
    ["not a url", DEV_SERVER, false],
  ])("navigating to %s stays in the application at %s: %s", (url, appUrl, allowed) => {
    const contents = fakeContents(appUrl);

    expect(contents.navigate(url)).toBe(!allowed);
  });

  test("a web view is never attached", () => {
    expect(fakeContents(PACKAGED).attachWebview()).toBe(true);
  });
});

describe("security senders", () => {
  test("the top-level frame of the application is trusted", () => {
    expect(isTrustedSender({ url: `${PACKAGED}runs/r1/log`, parent: null }, PACKAGED)).toBe(true);
  });

  test("a frame of another app scheme host is not trusted", () => {
    expect(isTrustedSender({ url: "app://other/runs/r1/log", parent: null }, PACKAGED)).toBe(false);
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

describe("security new windows", () => {
  test.each([
    ["https://doi.org/10.1093/ve/vex042", true],
    ["http://example.org", true],
    ["file:///etc/passwd", false],
    ["javascript:alert(1)", false],
  ])("a new window for %s is denied and opens in the system browser: %s", (url, external) => {
    const contents = fakeContents(PACKAGED);

    expect(contents.openWindow(url)).toStrictEqual({ action: "deny", opened: external ? [url] : [] });
  });
});

interface Preventable {
  preventDefault(): void;
}

function fakeContents(appUrl: string) {
  let openHandler: ((details: { url: string }) => { action: "deny" }) | undefined;
  const listeners = new Map<string, (event: Preventable, url: string) => void>();
  const opened: string[] = [];

  confineNavigation(
    {
      setWindowOpenHandler(handler) {
        openHandler = handler;
      },
      on(event: string, listener: (event: Preventable, url: string) => void) {
        listeners.set(event, listener);
      },
    },
    appUrl,
    (url) => {
      opened.push(url);
    },
  );

  const dispatch = (event: string, url: string): boolean => {
    let prevented = false;
    listeners.get(event)?.(
      {
        preventDefault: () => {
          prevented = true;
        },
      },
      url,
    );

    return prevented;
  };

  return {
    navigate: (url: string) => dispatch("will-navigate", url),
    attachWebview: () => dispatch("will-attach-webview", ""),
    openWindow: (url: string) => ({ action: openHandler?.({ url }).action, opened }),
  };
}
