import { baseLocale, type Locale } from "@neherlab/app-ui/i18n";

import { docsPath } from "./urls";
import type { PageManifest } from "./versions";

export function versionPageUrl(
  version: string,
  manifest: PageManifest | undefined,
  locale: Locale,
  page: string,
): string {
  const listed = manifest?.pages[locale]?.includes(page) === true;

  return docsPath(`/${version}`, listed ? page.split("/") : [], listed ? locale : baseLocale);
}
