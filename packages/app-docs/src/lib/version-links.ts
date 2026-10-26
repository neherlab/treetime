import { baseLocale, type Locale } from "@neherlab/app-ui/i18n";

import { docsPath } from "./urls";
import type { PageManifest } from "./versions";

export function versionPageUrl(
  version: string,
  manifest: PageManifest | undefined,
  locale: Locale,
  page: string,
): string {
  const pages = manifest?.pages[locale];

  if (pages?.includes(page) === true) {
    return docsPath(`/${version}`, page.split("/"), locale);
  }

  return docsPath(`/${version}`, [], pages?.includes("") === true ? locale : baseLocale);
}
