import { ThemeCycleButton } from "@neherlab/app-ui/chrome/ThemeCycleButton";
import { LocaleProvider, type Locale } from "@neherlab/app-ui/i18n";
import { ThemeProvider } from "@neherlab/app-ui/theme";

import { DocsSearch } from "./DocsSearch";
import { VersionSelect } from "./VersionSelect";

export function HeaderActions({
  locale,
  version,
  versionLabel,
  page,
}: {
  locale: Locale;
  version: string;
  versionLabel: string;
  page: string;
}) {
  return (
    <LocaleProvider locale={locale}>
      <ThemeProvider>
        <div className="ml-auto flex min-w-0 items-center gap-1 sm:gap-2">
          <DocsSearch />
          <VersionSelect version={version} label={versionLabel} page={page} />
          <ThemeCycleButton />
        </div>
      </ThemeProvider>
    </LocaleProvider>
  );
}
