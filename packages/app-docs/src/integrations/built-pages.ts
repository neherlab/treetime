import { locales, type Locale } from "@neherlab/app-ui/i18n";

import { assetPath } from "../lib/urls";

const NOT_FOUND_PAGE = "404";

export function pageUrl(base: string, file: string): string {
  return assetPath(base, withoutHtml(file));
}

export function isContentPage(file: string): boolean {
  return pageEntry(file) !== undefined;
}

export function manifestPages(files: readonly string[]): Record<string, string[]> {
  const entries = files.flatMap((file) => pageEntry(file) ?? []);

  return Object.fromEntries(
    locales.map((locale) => [
      locale,
      entries.flatMap((entry) => (entry.locale === locale ? [entry.page] : [])).toSorted(),
    ]),
  );
}

function pageEntry(file: string): PageEntry | undefined {
  if (!file.endsWith(".html")) {
    return undefined;
  }

  const [first = "", ...rest] = withoutHtml(file).split("/");
  const locale = locales.find((candidate) => candidate === first);
  const page = rest.join("/");

  return locale === undefined || page === NOT_FOUND_PAGE ? undefined : { locale, page };
}

interface PageEntry {
  locale: Locale;
  page: string;
}

function withoutHtml(file: string): string {
  return file.replace(/\.html$/u, "");
}
