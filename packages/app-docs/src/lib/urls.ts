import { baseLocale, locales, type Locale } from "@neherlab/app-ui/i18n";

export function docsUrl(slugs: readonly string[], locale: Locale): string {
  return docsPath(import.meta.env.BASE_URL, slugs, locale);
}

export function assetUrl(path: string): string {
  return assetPath(import.meta.env.BASE_URL, path);
}

export function docsPath(base: string, slugs: readonly string[], locale: Locale): string {
  const segments = [locale, ...slugs].filter((segment) => segment !== "");

  return `${trimTrailingSlashes(base)}/${segments.join("/")}`;
}

export function assetPath(base: string, path: string): string {
  return `${trimTrailingSlashes(base)}/${path.replace(/^\/+/u, "")}`;
}

export function docsLocale(value: string | undefined): Locale {
  return locales.find((locale) => locale === value) ?? baseLocale;
}

function trimTrailingSlashes(path: string): string {
  return path.replace(/\/+$/u, "");
}
