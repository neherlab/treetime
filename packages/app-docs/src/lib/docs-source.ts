import type { AppCommand } from "@neherlab/app-contracts";
import { baseLocale, locales, type Locale } from "@neherlab/app-ui/i18n";
import { defineI18n } from "fumadocs-core/i18n";
import { loader, type StaticSource } from "fumadocs-core/source";

import { docsLocale, docsPath } from "./urls";

export const docsI18n = defineI18n({
  languages: [...locales],
  defaultLanguage: baseLocale,
  parser: "dir",
  fallbackLanguage: baseLocale,
});

export function createDocsSource(files: DocsFiles, base: string) {
  const source: StaticSource<{ pageData: DocsPageData; metaData: DocsMetaData }> = {
    files: [
      ...files.pages.map((page) => ({ type: "page" as const, path: page.path, data: page.data })),
      ...files.metas.map((meta) => ({ type: "meta" as const, path: meta.path, data: meta.data })),
    ],
  };

  return loader({
    source,
    baseUrl: "/",
    i18n: docsI18n,
    url: (slugs, locale) => docsPath(base, slugs, docsLocale(locale)),
  });
}

export type DocsSource = ReturnType<typeof createDocsSource>;

export function isFallbackPage(page: { path: string }, locale: Locale): boolean {
  return locale !== baseLocale && !page.path.startsWith(`${locale}/`);
}

export interface DocsFiles {
  pages: readonly DocsPageFile[];
  metas: ReadonlyArray<{ path: string; data: DocsMetaData }>;
}

export interface DocsPageFile {
  path: string;
  data: DocsPageData;
}

export interface DocsPageData {
  title: string;
  description?: string | undefined;
  template: "doc" | "splash";
  tableOfContents: boolean;
  body: PageBody;
  editPath: string | undefined;
}

export type PageBody = { kind: "entry"; key: string } | { kind: "settings"; command: AppCommand };

export interface DocsMetaData {
  title?: string | undefined;
  pages?: string[] | undefined;
}
