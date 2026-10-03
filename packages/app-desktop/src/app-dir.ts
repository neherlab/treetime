import * as path from "node:path";

export const APP_DIR_ENV = "TREETIME_APP_DIR";

const CHECKOUT_APP_DIRS = {
  dev: { env: "TREETIME_DESKTOP_DEV_DIR", fallback: "tmp/app/treetime-dev" },
  prod: { env: "TREETIME_DESKTOP_PROD_DIR", fallback: "tmp/app/treetime-prod" },
} as const;

type CheckoutMode = keyof typeof CHECKOUT_APP_DIRS;

export function checkoutAppDir(
  env: Readonly<Record<string, string | undefined>>,
  mode: CheckoutMode,
  cwd: string,
): string {
  const { env: modeEnv, fallback } = CHECKOUT_APP_DIRS[mode];

  return path.resolve(cwd, nonEmpty(env[APP_DIR_ENV]) ?? nonEmpty(env[modeEnv]) ?? fallback);
}

function nonEmpty(value: string | undefined): string | undefined {
  return value === "" ? undefined : value;
}
