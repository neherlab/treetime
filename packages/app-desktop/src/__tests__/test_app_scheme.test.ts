import { describe, expect, test } from "vitest";

import { APP_URL, resolveAppAsset } from "../app-scheme";

const ROOT = "/opt/treetime/resources/app.asar/dist";

const INDEX = { kind: "file", path: `${ROOT}/index.html` };

const NOT_FOUND = { kind: "not-found" };

describe("app scheme", () => {
  test("the application URL opens the root route", () => {
    expect(new URL(APP_URL).pathname).toBe("/");
  });

  test.each([
    ["root", APP_URL],
    ["route", "app://treetime/new"],
    ["nested route", "app://treetime/runs/0123abcd/results"],
    ["route with query", "app://treetime/compare/a/b?tab=log"],
  ])("a %s serves the page", (_name, url) => {
    expect(resolveAppAsset(url, ROOT)).toStrictEqual(INDEX);
  });

  test.each([
    ["script", "app://treetime/assets/index-B1.js", "assets/index-B1.js"],
    ["style", "app://treetime/assets/index-C2.css", "assets/index-C2.css"],
    ["encoded name", "app://treetime/fonts/a%20b.woff2", "fonts/a b.woff2"],
  ])("a %s serves its file", (_name, url, file) => {
    expect(resolveAppAsset(url, ROOT)).toStrictEqual({ kind: "file", path: `${ROOT}/${file}` });
  });

  test.each([
    ["encoded parent directories", "app://treetime/..%2f..%2fsecret.txt"],
    ["encoded backslash parents", "app://treetime/..%5c..%5csecret.txt"],
    ["another host", "app://other/assets/index-B1.js"],
    ["another scheme", "file:///opt/treetime/resources/app.asar/dist/index.html"],
    ["malformed encoding", "app://treetime/assets/%E0.js"],
  ])("%s is not served", (_name, url) => {
    expect(resolveAppAsset(url, ROOT)).toStrictEqual(NOT_FOUND);
  });
});
