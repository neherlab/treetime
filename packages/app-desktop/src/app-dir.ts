import * as path from "node:path";

const CHECKOUT_EXAMPLES_DIR = "data";

const CHECKOUT_APP_DIRS = {
  dev: { env: "TREETIME_DESKTOP_DEV_DIR", fallback: "tmp/app/treetime-dev" },
  prod: { env: "TREETIME_DESKTOP_PROD_DIR", fallback: "tmp/app/treetime-prod" },
} as const;

type CheckoutMode = keyof typeof CHECKOUT_APP_DIRS;

type Env = Readonly<Record<string, string | undefined>>;

export interface AppFolderLaunch {
  appDirSwitch: string;
  env: Env;
  launchDir: string;
  checkout: { mode: CheckoutMode; root: string } | undefined;
}

export function appFolderEnv({ appDirSwitch, env, launchDir, checkout }: AppFolderLaunch): Record<string, string> {
  const appDir = nonEmpty(appDirSwitch);
  const fromSwitch = appDir === undefined ? {} : { TREETIME_APP_DIR: path.resolve(launchDir, appDir) };

  return checkout === undefined ? fromSwitch : checkoutEnv({ ...env, ...fromSwitch }, checkout.mode, checkout.root);
}

export function checkoutEnv(env: Env, mode: CheckoutMode, cwd: string) {
  const { env: modeEnv, fallback } = CHECKOUT_APP_DIRS[mode];

  return {
    TREETIME_APP_DIR: path.resolve(cwd, nonEmpty(env["TREETIME_APP_DIR"]) ?? nonEmpty(env[modeEnv]) ?? fallback),
    TREETIME_EXAMPLES_DIR: path.resolve(cwd, nonEmpty(env["TREETIME_EXAMPLES_DIR"]) ?? CHECKOUT_EXAMPLES_DIR),
  };
}

function nonEmpty(value: string | undefined): string | undefined {
  return value === "" ? undefined : value;
}
