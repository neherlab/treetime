import { createHash } from "node:crypto";
import { readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";

import mdx from "@astrojs/mdx";
import react from "@astrojs/react";
import sitemap from "@astrojs/sitemap";
import { settingCatalog } from "@neherlab/app-contracts/catalog";
import { contentSecurityPolicyDirectives } from "@neherlab/app-ui/build/content-security-policy-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import { baseLocale, locales } from "@neherlab/app-ui/i18n";
import tailwindcss from "@tailwindcss/vite";
import expressiveCode from "astro-expressive-code";
import { defineConfig } from "astro/config";

import { postBuild } from "./src/integrations/post-build";
import { versionListSchema } from "./src/lib/versions";
import { docsMarkdownProcessor } from "./src/markdown/pipeline";
import { THEME_SCRIPT } from "./src/theme/theme-script";
import versionsFile from "./versions.json";

const DOCS_VERSION = process.env["TREETIME_DOCS_VERSION"] ?? "dev";

const SITE = process.env["TREETIME_DOCS_SITE_URL"] ?? "https://docs.treetime.neherlab.click";

const BASE = `/${DOCS_VERSION}`;

const VERSIONS = versionListSchema.parse(versionsFile);

const SETTINGS_PAGES = new URL("src/content/docs/reference/settings/", import.meta.url);

const HASHED_ASSET = "assets/[name]-[hash][extname]";

const HASHED_SCRIPT = "assets/[name]-[hash].js";

const DIRECTIVES = contentSecurityPolicyDirectives(false, { webAssembly: true });

checkVersion();

checkSettingsPages();

export default defineConfig({
  site: SITE,
  base: BASE,
  trailingSlash: "never",
  publicDir,
  build: { format: "file", assets: "assets" },
  i18n: {
    locales: [...locales],
    defaultLocale: baseLocale,
    routing: { prefixDefaultLocale: false },
  },
  server: ({ command }) => ({
    port: Number(
      command === "preview"
        ? (process.env["TREETIME_DOCS_SERVE_PORT"] ?? "6700")
        : (process.env["TREETIME_DOCS_PORT"] ?? "6100"),
    ),
    allowedHosts: [".localhost"],
  }),
  security: {
    csp: {
      directives: DIRECTIVES.filter(isOtherDirective).map(directiveText),
      scriptDirective: {
        resources: [...directiveSources("script-src")],
        hashes: [`sha256-${createHash("sha256").update(THEME_SCRIPT).digest("base64")}`],
      },
      styleDirective: { resources: [...directiveSources("style-src")] },
    },
  },
  markdown: {
    syntaxHighlight: false,
    processor: docsMarkdownProcessor(BASE),
  },
  integrations: [
    expressiveCode({
      themes: ["everforest-light", "everforest-dark"],
      useDarkModeMediaQuery: false,
      themeCssSelector: (theme) => (theme.type === "dark" ? ".dark" : false),
      styleOverrides: {
        borderRadius: "var(--radius)",
        borderColor: "var(--border)",
        codeFontFamily: '"IBM Plex Mono", ui-monospace, monospace',
        uiFontFamily: '"Lato", ui-sans-serif, sans-serif',
        frames: { frameBoxShadowCssValue: "none" },
      },
    }),
    mdx({ extendMarkdownConfig: true }),
    react(),
    sitemap({
      i18n: { defaultLocale: baseLocale, locales: Object.fromEntries(locales.map((locale) => [locale, locale])) },
    }),
    postBuild({ base: BASE, version: DOCS_VERSION, versions: VERSIONS, site: SITE }),
  ],
  vite: {
    plugins: [icons(), tailwindcss()],
    build: { rolldownOptions: { output: { assetFileNames: HASHED_ASSET } } },
    environments: {
      client: {
        build: {
          rolldownOptions: {
            output: { assetFileNames: HASHED_ASSET, chunkFileNames: HASHED_SCRIPT, entryFileNames: HASHED_SCRIPT },
          },
        },
      },
    },
  },
});

function isOtherDirective(directive: (typeof DIRECTIVES)[number]): directive is OtherDirective {
  return directive[0] !== "script-src" && directive[0] !== "style-src";
}

function directiveText([name, sources]: OtherDirective): `${OtherDirective[0]} ${string}` {
  return `${name} ${sources.join(" ")}`;
}

function directiveSources(name: DirectiveName): readonly string[] {
  return DIRECTIVES.find(([directive]) => directive === name)?.[1] ?? [];
}

function checkVersion(): void {
  if (!VERSIONS.some((entry) => entry.version === DOCS_VERSION)) {
    throw new Error(
      `The docs version ${DOCS_VERSION} has no entry in versions.json, which lists ${VERSIONS.map((entry) => entry.version).join(", ")}`,
    );
  }
}

function checkSettingsPages(): void {
  const pages = new Set(
    readdirSync(fileURLToPath(SETTINGS_PAGES))
      .filter((file) => file.endsWith(".mdx"))
      .map((file) => file.replace(/\.mdx$/u, "")),
  );

  const commands = new Set<string>(settingCatalog.commands.map((command) => command.command));
  const missing = [...commands].filter((command) => !pages.has(command));
  const extra = [...pages].filter((page) => !commands.has(page));

  if (missing.length > 0 || extra.length > 0) {
    throw new Error(
      `The settings reference pages in ${fileURLToPath(SETTINGS_PAGES)} do not match the setting catalog: ` +
        `missing pages for ${missing.join(", ") || "none"}, pages without a command: ${extra.join(", ") || "none"}`,
    );
  }
}

type DirectiveName = (typeof DIRECTIVES)[number][0];

type OtherDirective = readonly [Exclude<DirectiveName, "script-src" | "style-src">, readonly string[]];
