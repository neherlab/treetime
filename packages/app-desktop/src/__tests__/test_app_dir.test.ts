import { describe, expect, test } from "vitest";

import { checkoutAppDir } from "../app-dir";

describe("app_dir checkout", () => {
  test.each([
    ["dev", "/checkout/tmp/app/treetime-dev"],
    ["prod", "/checkout/tmp/app/treetime-prod"],
  ] as const)("the %s mode keeps its data in its own folder of the checkout", (mode, expected) => {
    expect(checkoutAppDir({}, mode, "/checkout")).toBe(expected);
  });

  test("the folder of the mode comes from its environment variable, relative to the checkout", () => {
    expect(checkoutAppDir({ TREETIME_DESKTOP_DEV_DIR: "tmp/app/mine" }, "dev", "/checkout")).toBe(
      "/checkout/tmp/app/mine",
    );
  });

  test("the variable of the other mode does not apply", () => {
    expect(checkoutAppDir({ TREETIME_DESKTOP_DEV_DIR: "/elsewhere" }, "prod", "/checkout")).toBe(
      "/checkout/tmp/app/treetime-prod",
    );
  });

  test("the shared TREETIME_APP_DIR applies to both modes and wins over the variable of the mode", () => {
    expect(
      checkoutAppDir(
        { TREETIME_APP_DIR: "/data/treetime", TREETIME_DESKTOP_PROD_DIR: "/elsewhere" },
        "prod",
        "/checkout",
      ),
    ).toBe("/data/treetime");
  });

  test("an empty variable counts as unset", () => {
    expect(checkoutAppDir({ TREETIME_APP_DIR: "", TREETIME_DESKTOP_DEV_DIR: "" }, "dev", "/checkout")).toBe(
      "/checkout/tmp/app/treetime-dev",
    );
  });
});
