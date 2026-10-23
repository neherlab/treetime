import { Glob, plugin, TOML } from "bun";
import * as fc from "fast-check";
import * as z from "zod";

import { AUSPICE_EXTENSIONS, AUSPICE_SOURCE, transformAuspiceSource } from "./packages/app-ui/build/auspice-vite.ts";

const ROOT = import.meta.dir;

const ICON_SOURCES = "packages/*/{src,renderer}/**/*.{ts,tsx}";

const ICON_IMPORT = /from "(?<specifier>~icons\/[^"]+)"/gu;

const bunfigSchema = z.object({ test: z.object({ seed: z.number() }) });

const { test: testConfig } = bunfigSchema.parse(TOML.parse(await Bun.file(`${ROOT}/bunfig.toml`).text()));

fc.configureGlobal({ ...fc.readConfigureGlobal(), seed: testConfig.seed });

const iconSpecifiers = await importedIcons();

const loadIcon = lazyIconLoader();

plugin({
  name: "treetime-test-modules",
  setup(build) {
    /* HACK: Bun runtime plugins never call onResolve for bare specifiers (oven-sh/bun#40579), so each icon
       import found in the sources becomes a virtual module instead of being resolved on demand. */
    for (const specifier of iconSpecifiers) {
      build.module(specifier, async () => ({ contents: await loadIcon(specifier), loader: "jsx" }));
    }

    build.module(AUSPICE_EXTENSIONS, () => ({ contents: "export default {};", loader: "js" }));

    build.onLoad({ filter: AUSPICE_SOURCE }, async ({ path }) => ({
      contents: (await transformAuspiceSource(await Bun.file(path).text(), path)).code,
      loader: "js",
    }));
  },
});

async function importedIcons(): Promise<Set<string>> {
  const specifiers = new Set<string>();

  for await (const file of new Glob(ICON_SOURCES).scan({ cwd: ROOT, absolute: true })) {
    for (const match of (await Bun.file(file).text()).matchAll(ICON_IMPORT)) {
      const specifier = match.groups?.["specifier"];

      if (specifier !== undefined) specifiers.add(specifier);
    }
  }

  return specifiers;
}

function lazyIconLoader(): (specifier: string) => Promise<string> {
  let loader: Promise<(specifier: string) => Promise<string>> | undefined;

  return async (specifier) => {
    loader ??= import("./packages/app-ui/build/icons-vite.ts").then(({ iconModuleLoader }) => iconModuleLoader());

    return (await loader)(specifier);
  };
}
