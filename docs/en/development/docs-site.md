---
title: Docs site
description: How the docs site is built, where its pages are, and how it shares the app interface.
---

The docs site (`packages/app-docs`) is a prerendered multi-page site in the look of the app. It is its own [Astro](https://astro.build) site: the headless library [`fumadocs-core`](https://www.fumadocs.dev) supplies the page tree, the language fallback, the neighbour links, and the tracking of the visible heading, and `@neherlab/app-ui` renders every visible part. The docs are a separate build, not a route of the app, because they change once per release while the app deploys on every commit.

| Task                                                                       | Command           |
| -------------------------------------------------------------------------- | ----------------- |
| Start the dev server for writing pages, with live reload                   | `just docs`       |
| Build the docs as published and serve the build like the docs host         | `just docs-prod`  |
| Build into `packages/app-docs/dist/` and `packages/app-docs/dist-root/`    | `just build-docs` |

`just docs` runs in the `sess` session `tt-dev-docs` and `just docs-prod` in `tt-dev-docs-prod`. The dev server serves `http://127.0.0.1:<port>/dev` on the docs port of the checkout (14005 in the main checkout), which leads to the English home page `/dev/en`, and reloads pages and styles as they change; it runs no post-build step, so search and the version list are missing there. `just docs-prod` builds the site and serves it on the docs serve port (14006 in the main checkout) with the file layout of the docs host: the site root files (`versions.json`, the page that leads to the preferred version) at `/`, the build at `/dev`, pages without their `.html` extension, and the not-found page of the locale for a missing page (`packages/app-docs/src/preview/serve-site.ts`). `TREETIME_DOCS_PORT` and `TREETIME_DOCS_SERVE_PORT` override the ports, and `just status` prints them.

The app and the docs share three layers, and each surface wires the rest itself:

- **Styles**: `packages/app-docs/src/styles/docs.css` loads Tailwind, `tokens.css`, `components.css`, and `base.css` of app-ui, and adds the prose styles of `@tailwindcss/typography` with the colors of the app tokens. An app-ui component inside an article sits in a `not-prose` element, so it looks as it does in the app
- **Components**: the header parts (`packages/app-ui/src/chrome/`), the theme control, the search palette, the sidebar items, the alerts, and the setting presentation (`packages/app-ui/src/settings/help.tsx`) are app-ui components. The lint bans `@base-ui/react`, `react-aria-components`, and `next-themes` in the docs, so a missing control goes into `packages/app-ui/src/ui/`
- **Messages**: the strings of the shared parts and of the docs layout are Paraglide messages of app-ui (see UI messages below)
- **Surface wiring**: the app keeps its router, API, and stores; the docs wire the content collections, `fumadocs-core`, Pagefind, and the version data

The layout (`src/layouts/DocsLayout.astro`) has the header (`src/components/DocsHeader.astro`), the navigation tree, the article, and the table of contents. Astro renders each React component as its own tree, and a `client:` directive works only on a component that the `.astro` or `.mdx` file imports itself, so every interactive part is one island: the header actions (search, version select, theme), the mobile navigation, the table of contents, the version banner, the settings reference, and the interactive page components. Everything else renders at build time and ships no JavaScript. `src/lib/urls.ts` builds every link of the site. Logic lives in `.ts` and `.tsx` files, and `.astro` files stay thin, because the type check does not read `.astro` files.

`astro build` alone produces the complete output: one post-build step (`src/integrations/post-build.ts`) writes the Pagefind search index under the canonical page URLs, the page manifest `pages.json`, the site root files in `dist-root/`, and the `.gz` copies, and checks that each executable inline script of each page has its hash in the content security policy.

## Pages

The pages live in `docs/` at the top of the repository, apart from the code of the site, so authors and translators work in one shallow folder. `src/content.config.ts` reads them as the content collection `docs`:

- **One folder per language**: `docs/en/` and `docs/de/`, with the same paths. A German path without a German file shows the English page with a notice
- **Files**: `.md`, or `.mdx` for a page with components. `meta.json` files set the title of a folder and the order of its pages
- **Callouts** are directives (`:::note`, `:::tip`, `:::warning`, `:::danger`); another type stops the build. A heading that starts with `1. ` becomes a step
- **Links** are root-relative without the version and the language (`/showcase/math`); the Markdown pipeline (`src/markdown/pipeline.ts`) adds both
- **Math** uses KaTeX with `\tag{n}` and an HTML anchor before the display block, because KaTeX has no `\label`; shared macros are in `src/math/macros.ts`
- **Edit links**: each page links to its file in `docs/` on GitHub

The showcase pages hold placeholder text that exercises each feature; real guides replace them. `docs/en/showcase/writing.mdx` lists the syntax of each feature.

## Components in pages

Pages use components without an import. `src/markdown/remark-provided-components.ts` adds the import of each component that a page uses, from the list in `src/markdown/provided-imports.ts`, and adds the `client:` directive of the interactive ones (`Disclosure`, `SyncedTabs`), because Astro hydrates only components that the page imports itself. An element named `Icon<Set><Name>`, such as `IconLucideClock`, becomes an import of `~icons/<set>/<name>` through the resolver of `unplugin-icons`. To provide a new component, add it to `ProvidedComponents` in `src/components/provided-components.ts` and to the list; the type check fails when the two differ.

Editors get the component names and props in `.mdx` files from the MDX language server (the MDX extension of VS Code): `docs/tsconfig.json` turns on its type check (`mdx.checkMdx`), loads the remark plugins that parse front matter, directives, and math, and includes the global type `MDXProvidedComponents`.

## Generated pages

Pages that come from the Rust code have no file in `docs/`. The settings pages render the setting catalog, one page per command, and `src/lib/source.ts` adds them to the page tree. The command-line reference is the collection `generated` in `packages/app-docs/src/generated/en/`, written by `just gen cli-docs`.

## URLs

Every page URL names the version and the language: `/<version>/<language>/<page>`, such as `/dev/en/showcase/math` and `/dev/de`. `/<version>` leads to the English home page, and the site root `/` to the English home page of the preferred version. `src/lib/urls.ts` builds every link of the site.

## Translations

Translation agents write the pages of `docs/de/` from the English pages. A translation records the hash of the English file it was translated from: `sourceHash` in the front matter of a page, or a `sourceHash` key in a `meta.json` file. `just docs-i18n-status` lists one line per translation to write, with the hash to record:

- `missing`: the English file has no translation
- `stale`: the English file changed after the translation, or the translation records no hash
- `orphan`: the translation has no English file; delete it or move it with its English page

A translation keeps the structure of its English page: the same headings, components, links, and code blocks, with only the text translated. The developer pages are translated like every other page.

## UI messages

The UI strings of the shared parts, the docs, and (in English) the app are Paraglide messages: `packages/app-ui/messages/en.json` and `de.json`, in the inlang project `packages/app-ui/project.inlang`, whose `settings.json` loads its plugins from `node_modules`. `just gen i18n` compiles them into `packages/app-ui/src/i18n/generated/`, which the builds, the tests, and the type check read; no build runs the Paraglide Vite plugin. To add a message, add the key to both files and run `just gen i18n`; a test fails when a key or a parameter is missing in a locale. A component reads its locale with `useLocale()` and passes it to each call, `m.<key>(params, { locale })`, so pages that render side by side at build time cannot mix languages. A new locale goes into `settings.json`, a message file, and the `i18n` and route settings of the docs. Setting texts, help texts, CLI help, errors, and run warnings come from Rust and stay English on every page.
