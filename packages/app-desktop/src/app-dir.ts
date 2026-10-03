import * as path from "node:path";

const CHECKOUT_EXAMPLES_DIR = "data";

const CHECKOUT_APP_DIRS = {
  dev: { env: "TREETIME_DESKTOP_DEV_DIR", fallback: "tmp/app/treetime-dev" },
  prod: { env: "TREETIME_DESKTOP_PROD_DIR", fallback: "tmp/app/treetime-prod" },
} as const;

type CheckoutMode = keyof typeof CHECKOUT_APP_DIRS;

export function checkoutEnv(env: Readonly<Record<string, string | undefined>>, mode: CheckoutMode, cwd: string) {
  const { env: modeEnv, fallback } = CHECKOUT_APP_DIRS[mode];

  return {
    TREETIME_APP_DIR: path.resolve(cwd, nonEmpty(env["TREETIME_APP_DIR"]) ?? nonEmpty(env[modeEnv]) ?? fallback),
    TREETIME_EXAMPLES_DIR: path.resolve(cwd, nonEmpty(env["TREETIME_EXAMPLES_DIR"]) ?? CHECKOUT_EXAMPLES_DIR),
  };
}

function nonEmpty(value: string | undefined): string | undefined {
  return value === "" ? undefined : value;
}
