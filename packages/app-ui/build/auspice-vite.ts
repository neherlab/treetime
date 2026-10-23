import { readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";

import type { Plugin } from "vite";

export const AUSPICE_SOURCE =
  /[\\/]node_modules[\\/](?:\.bun[\\/][^\\/]+[\\/]node_modules[\\/])?auspice[\\/]src[\\/][^?]+\.js(?:\?.*)?$/u;

const APP_UI_PACKAGE = "@neherlab/app-ui";

export const AUSPICE_EXTENSIONS = "@extensions";

const AUSPICE_EXTENSIONS_MODULE = fileURLToPath(new URL("../src/auspice/extensions.ts", import.meta.url));

const AUSPICE_DECLARATIONS = fileURLToPath(new URL("../src/auspice/types", import.meta.url));

const DECLARATION = /\.d\.ts$/u;

const INDEX_MODULE = /\/index$/u;

const SHARED_PACKAGES = ["react", "react-dom"];

const AUSPICE_TIMERS = /[\\/]auspice[\\/]src[\\/]util[\\/]perf\.js(?:\?.*)?$/u;

const AUSPICE_TIMERS_STUB = "export const timerStart = () => {};\nexport const timerEnd = () => {};\n";

const AUSPICE_TRANSFORM = {
  filter: { id: AUSPICE_SOURCE },
  handler: transformAuspiceSource,
} satisfies Plugin["transform"];

export async function transformAuspiceSource(code: string, id: string) {
  if (AUSPICE_TIMERS.test(id)) {
    return { code: AUSPICE_TIMERS_STUB, map: null };
  }

  const { transformWithOxc } = await import("vite");

  const result = await transformWithOxc(code, id, {
    lang: "jsx",
    jsx: { runtime: "automatic" },
    decorator: { legacy: true },
    sourcemap: true,
  });

  return { code: result.code, map: result.map ?? null };
}

export function auspice(): Plugin {
  return {
    name: "treetime-auspice",
    enforce: "pre",
    config() {
      return {
        define: {
          "process.env.SKIP_REDUX_CHECKS": JSON.stringify(""),
          "process.env.ENABLE_SERVICE_WORKER": "false",
        },
        resolve: { alias: { [AUSPICE_EXTENSIONS]: AUSPICE_EXTENSIONS_MODULE }, dedupe: SHARED_PACKAGES },
        optimizeDeps: {
          extensions: [".tsx"],
          include: auspiceEntries(readdirSync(AUSPICE_DECLARATIONS, { recursive: true, encoding: "utf8" })).map(
            (entry) => `${APP_UI_PACKAGE} > ${entry}`,
          ),
          rolldownOptions: { plugins: [{ name: "treetime-auspice-source", transform: AUSPICE_TRANSFORM }] },
        },
      };
    },
    transform: AUSPICE_TRANSFORM,
  };
}

export function auspiceEntries(declarations: readonly string[]): string[] {
  return declarations
    .map((path) => path.replaceAll("\\", "/"))
    .filter((path) => DECLARATION.test(path))
    .map((path) => `auspice/src/${path.replace(DECLARATION, "").replace(INDEX_MODULE, "")}`)
    .toSorted();
}
