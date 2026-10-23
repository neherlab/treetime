import { baseLocale, locales, type Locale } from "@neherlab/app-ui/i18n";
import { defineI18n } from "fumadocs-core/i18n";
import { loader, type StaticSource } from "fumadocs-core/source";

import { docsLocale, docsPath } from "./urls";

export const CONTENT_ROOT = "src/content/docs/";

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
  pages: ReadonlyArray<{ path: string; data: DocsPageData }>;
  metas: ReadonlyArray<{ path: string; data: DocsMetaData }>;
}

export interface DocsPageData {
  title: string;
  description?: string | undefined;
  template: "doc" | "splash";
  tableOfContents: boolean;
  editUrl?: false | undefined;
  entryId: string;
  file: string;
}

export interface DocsMetaData {
  title?: string | undefined;
  pages?: string[] | undefined;
}
