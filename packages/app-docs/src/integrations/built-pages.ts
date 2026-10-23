import { baseLocale, locales, type Locale } from "@neherlab/app-ui/i18n";

import { assetPath, docsPath } from "../lib/urls";

const NOT_FOUND_PAGE = "404";

export function pageUrl(base: string, file: string): string {
  const path = withoutHtml(file);

  return path === "index" ? docsPath(base, [], baseLocale) : assetPath(base, path);
}

export function isContentPage(file: string): boolean {
  const path = withoutHtml(file);

  return file.endsWith(".html") && path !== NOT_FOUND_PAGE && !path.endsWith(`/${NOT_FOUND_PAGE}`);
}

export function manifestPages(files: readonly string[]): Record<string, string[]> {
  const entries = files.flatMap((file) => (isContentPage(file) ? [pageEntry(withoutHtml(file))] : []));

  return Object.fromEntries(
    locales.map((locale) => [
      locale,
      entries.flatMap((entry) => (entry.locale === locale ? [entry.page] : [])).toSorted(),
    ]),
  );
}

function pageEntry(path: string): PageEntry {
  const locale =
    locales.find((candidate) => candidate !== baseLocale && (path === candidate || path.startsWith(`${candidate}/`))) ??
    baseLocale;

  const page = locale === baseLocale ? path : path.slice(locale.length + 1);

  return { locale, page: page === "index" ? "" : page };
}

interface PageEntry {
  locale: Locale;
  page: string;
}

function withoutHtml(file: string): string {
  return file.replace(/\.html$/u, "");
}
