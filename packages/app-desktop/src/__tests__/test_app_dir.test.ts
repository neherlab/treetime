import { describe, expect, test } from "vitest";

import { appFolderEnv, checkoutEnv } from "../app-dir";

describe("app_dir checkout", () => {
  test.each([
    ["dev", "/checkout/tmp/app/treetime-dev"],
    ["prod", "/checkout/tmp/app/treetime-prod"],
  ] as const)("the %s mode keeps its data in its own folder of the checkout", (mode, expected) => {
    expect(checkoutEnv({}, mode, "/checkout")).toStrictEqual({
      TREETIME_APP_DIR: expected,
      TREETIME_EXAMPLES_DIR: "/checkout/data",
    });
  });

  test("the folder of the mode comes from its environment variable, relative to the checkout", () => {
    expect(checkoutEnv({ TREETIME_DESKTOP_DEV_DIR: "tmp/app/mine" }, "dev", "/checkout").TREETIME_APP_DIR).toBe(
      "/checkout/tmp/app/mine",
    );
  });

  test("the variable of the other mode does not apply", () => {
    expect(checkoutEnv({ TREETIME_DESKTOP_DEV_DIR: "/elsewhere" }, "prod", "/checkout").TREETIME_APP_DIR).toBe(
      "/checkout/tmp/app/treetime-prod",
    );
  });

  test("the shared TREETIME_APP_DIR applies to both modes and wins over the variable of the mode", () => {
    const env = { TREETIME_APP_DIR: "/data/treetime", TREETIME_DESKTOP_PROD_DIR: "/elsewhere" };

    expect(checkoutEnv(env, "prod", "/checkout").TREETIME_APP_DIR).toBe("/data/treetime");
  });

  test("the examples come from the checkout data unless a variable names another folder", () => {
    expect(checkoutEnv({ TREETIME_EXAMPLES_DIR: "examples" }, "dev", "/checkout").TREETIME_EXAMPLES_DIR).toBe(
      "/checkout/examples",
    );
  });

  test("an empty variable counts as unset", () => {
    const env = { TREETIME_APP_DIR: "", TREETIME_DESKTOP_DEV_DIR: "", TREETIME_EXAMPLES_DIR: "" };

    expect(checkoutEnv(env, "dev", "/checkout")).toStrictEqual({
      TREETIME_APP_DIR: "/checkout/tmp/app/treetime-dev",
      TREETIME_EXAMPLES_DIR: "/checkout/data",
    });
  });
});

describe("app_dir switch", () => {
  test("the switch names the app folder relative to the launch folder of an installed app", () => {
    expect(
      appFolderEnv({ appDirSwitch: "check", env: {}, launchDir: "/home/user", checkout: undefined }),
    ).toStrictEqual({ TREETIME_APP_DIR: "/home/user/check" });
  });

  test("an installed app without the switch keeps the environment", () => {
    expect(
      appFolderEnv({
        appDirSwitch: "",
        env: { TREETIME_APP_DIR: "/data" },
        launchDir: "/home/user",
        checkout: undefined,
      }),
    ).toStrictEqual({});
  });

  test("the switch wins over TREETIME_APP_DIR and over the folder of the checkout mode", () => {
    const env = { TREETIME_APP_DIR: "/data/treetime", TREETIME_DESKTOP_DEV_DIR: "/elsewhere" };

    expect(
      appFolderEnv({
        appDirSwitch: "/tmp/check",
        env,
        launchDir: "/home/user",
        checkout: { mode: "dev", root: "/checkout" },
      }),
    ).toStrictEqual({ TREETIME_APP_DIR: "/tmp/check", TREETIME_EXAMPLES_DIR: "/checkout/data" });
  });

  test("a relative switch resolves from the launch folder, not from the checkout", () => {
    expect(
      appFolderEnv({
        appDirSwitch: "check",
        env: {},
        launchDir: "/home/user",
        checkout: { mode: "prod", root: "/checkout" },
      })["TREETIME_APP_DIR"],
    ).toBe("/home/user/check");
  });

  test("a checkout without the switch takes the folder of its mode", () => {
    expect(
      appFolderEnv({
        appDirSwitch: "",
        env: {},
        launchDir: "/home/user",
        checkout: { mode: "dev", root: "/checkout" },
      }),
    ).toStrictEqual({ TREETIME_APP_DIR: "/checkout/tmp/app/treetime-dev", TREETIME_EXAMPLES_DIR: "/checkout/data" });
  });
});
