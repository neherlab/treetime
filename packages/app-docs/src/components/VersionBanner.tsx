import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { Alert, AlertDescription } from "@neherlab/app-ui/ui/alert";
import History from "~icons/lucide/history";

import { versionPageUrl } from "../lib/version-links";
import { useManifest, useSiteVersions } from "../lib/version-state";
import { preferredVersion } from "../lib/versions";

export function VersionBanner({ version, page, locale }: { version: string; page: string; locale: Locale }) {
  const versions = useSiteVersions();
  const preferred = versions === undefined ? undefined : preferredVersion(versions);
  const outdated = preferred !== undefined && preferred.version !== version;
  const manifest = useManifest(outdated ? preferred.version : undefined);
  const current = versions?.find((entry) => entry.version === version);

  if (!outdated || current === undefined) {
    return null;
  }

  return (
    <LocaleProvider locale={locale}>
      <Alert role="note" variant="warning" className="not-prose mb-6" data-pagefind-ignore>
        <History aria-hidden />
        <AlertDescription>
          {m.docs_version_banner({ current: current.label, preferred: preferred.label }, { locale })}{" "}
          <a href={versionPageUrl(preferred.version, manifest, locale, page)}>
            {m.docs_version_banner_link({ preferred: preferred.label }, { locale })}
          </a>
        </AlertDescription>
      </Alert>
    </LocaleProvider>
  );
}
