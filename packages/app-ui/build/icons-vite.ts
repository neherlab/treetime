import { fileURLToPath } from "node:url";

import iconsUnplugin from "unplugin-icons";
import type { Options } from "unplugin-icons";
import Icons from "unplugin-icons/vite";
import type { PluginOption } from "vite";
import * as z from "zod";

const ICON_SIZE = "24";

const ICON_OPTIONS: Options = {
  compiler: "jsx",
  jsx: "react",
  collectionsNodeResolvePath: fileURLToPath(new URL("..", import.meta.url)),
  iconCustomizer(_collection, _icon, props) {
    props["width"] = ICON_SIZE;
    props["height"] = ICON_SIZE;
  },
};

const HOOK_CONTEXT: HookContext = {
  addWatchFile: () => undefined,
  emitFile: () => undefined,
  getWatchFiles: () => [],
  parse: () => {
    throw new Error("the icon hooks parse no modules outside a bundler");
  },
  error: (message: Parameters<HookContext["error"]>[0]) => {
    throw new Error("unplugin-icons reported an error", { cause: message });
  },
  warn: () => undefined,
};

const resolvedIconSchema = z.string();

const loadedIconSchema = z.union([z.string(), z.object({ code: z.string() }).transform(({ code }) => code)]);

export function icons(): PluginOption {
  return Icons(ICON_OPTIONS);
}

export function iconModuleLoader(): (specifier: string) => Promise<string> {
  const hooks = iconHooks();
  const resolveId = hookHandler(hooks.resolveId);
  const load = hookHandler(hooks.load);

  return async (specifier) => {
    const id = resolvedIconSchema.safeParse(
      await resolveId.call(HOOK_CONTEXT, specifier, undefined, { isEntry: false }),
    );

    if (!id.success) {
      throw new Error(`unplugin-icons does not resolve the icon import ${specifier}: ${z.prettifyError(id.error)}`);
    }

    const code = loadedIconSchema.safeParse(await load.call(HOOK_CONTEXT, id.data));

    if (!code.success) {
      throw new Error(`unplugin-icons has no icon for the import ${specifier}: ${z.prettifyError(code.error)}`);
    }

    return code.data;
  };
}

type IconHooks = Exclude<ReturnType<typeof iconsUnplugin.raw>, readonly unknown[]>;

type ResolveIdHook = Exclude<NonNullable<IconHooks["resolveId"]>, { readonly handler: CallableFunction }>;

type HookContext = ThisParameterType<ResolveIdHook>;

function iconHooks(): IconHooks {
  const hooks = iconsUnplugin.raw(ICON_OPTIONS, { framework: "bun", versions: {} });

  if (Array.isArray(hooks)) {
    throw new TypeError("unplugin-icons returned several plugins; the icon module loader expects one");
  }

  return hooks;
}

function hookHandler<T extends CallableFunction>(hook: T | { readonly handler: T } | undefined): T {
  if (hook === undefined) {
    throw new Error("unplugin-icons lacks the resolveId and load hooks that the icon module loader calls");
  }

  return "handler" in hook ? hook.handler : hook;
}
