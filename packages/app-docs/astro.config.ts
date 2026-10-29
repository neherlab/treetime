import { createHash } from "node:crypto";
import { fileURLToPath } from "node:url";

import mdx from "@astrojs/mdx";
import react from "@astrojs/react";
import sitemap from "@astrojs/sitemap";
import { contentSecurityPolicyDirectives } from "@neherlab/app-ui/build/content-security-policy-vite";
import { icons } from "@neherlab/app-ui/build/icons-vite";
import { publicDir } from "@neherlab/app-ui/build/public-dir";
import { baseLocale, locales } from "@neherlab/app-ui/i18n";
import tailwindcss from "@tailwindcss/vite";
import expressiveCode from "astro-expressive-code";
import { defineConfig } from "astro/config";

import { postBuild } from "./src/integrations/post-build";
import { CONTENT_ROOTS } from "./src/lib/content-roots";
import { versionListSchema } from "./src/lib/versions";
import { docsMarkdownProcessor } from "./src/markdown/pipeline";
import { THEME_SCRIPT } from "./src/theme/theme-script";
import versionsFile from "./versions.json";

const DOCS_VERSION = process.env["TREETIME_DOCS_VERSION"] ?? "dev";

const SITE = process.env["TREETIME_DOCS_SITE_URL"] ?? "https://docs.treetime.neherlab.click";

const BASE = `/${DOCS_VERSION}`;

const VERSIONS = versionListSchema.parse(versionsFile);

const CONTENT_ROOT_PATHS = Object.values(CONTENT_ROOTS).map((root) => fileURLToPath(new URL(root, import.meta.url)));

const HASHED_ASSET = "assets/[name]-[hash][extname]";

const HASHED_SCRIPT = "assets/[name]-[hash].js";

const DIRECTIVES = contentSecurityPolicyDirectives(false, { webAssembly: true });

checkVersion();

export default defineConfig({
  site: SITE,
  base: BASE,
  trailingSlash: "never",
  experimental: { contentIntellisense: true },
  devToolbar: { enabled: false },
  publicDir,
  build: { format: "file", assets: "assets" },
  i18n: {
    locales: [...locales],
    defaultLocale: baseLocale,
    routing: { prefixDefaultLocale: true, redirectToDefaultLocale: false },
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
    processor: docsMarkdownProcessor(BASE, CONTENT_ROOT_PATHS),
  },
  integrations: [
    expressiveCode({
      themes: ["github-light", "github-dark"],
      useDarkModeMediaQuery: false,
      themeCssSelector: (theme) => (theme.type === "dark" ? ".dark" : false),
      styleOverrides: {
        borderRadius: "var(--corner-lg)",
        borderColor: "var(--border)",
        codeBackground: "var(--muted)",
        codeFontFamily: "var(--typeface-mono)",
        uiFontFamily: "var(--typeface-sans)",
        frames: {
          frameBoxShadowCssValue: "none",
          editorActiveTabIndicatorTopColor: "var(--primary)",
          editorActiveTabBackground: "var(--muted)",
          editorTabBarBackground: "var(--card)",
          terminalTitlebarBackground: "var(--card)",
          terminalBackground: "var(--muted)",
        },
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

type DirectiveName = (typeof DIRECTIVES)[number][0];

type OtherDirective = readonly [Exclude<DirectiveName, "script-src" | "style-src">, readonly string[]];
