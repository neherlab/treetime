import { m, useLocale, type Locale } from "@neherlab/app-ui/i18n";
import { NativeSelect, NativeSelectOption } from "@neherlab/app-ui/ui/native-select";
import { useCallback, type ChangeEvent } from "react";

import { versionPageUrl } from "../lib/version-links";
import { manifestOf, useSiteVersions } from "../lib/version-state";

export function VersionSelect({ version, label, page }: { version: string; label: string; page: string }) {
  const locale = useLocale();
  const versions = useSiteVersions();

  const choose = useCallback(
    async (event: ChangeEvent<HTMLSelectElement>) => {
      const target = event.currentTarget.value;
      window.location.assign(await targetUrl(target, locale, page));
    },
    [locale, page],
  );

  const onChange = useCallback((event: ChangeEvent<HTMLSelectElement>) => void choose(event), [choose]);

  if (versions === undefined || versions.length < 2) {
    return <span className="text-muted-foreground hidden font-mono text-xs sm:inline">{label}</span>;
  }

  return (
    <NativeSelect
      size="sm"
      value={version}
      onChange={onChange}
      aria-label={m.docs_version({}, { locale })}
      className="hidden sm:block"
    >
      {versions.map((entry) => (
        <NativeSelectOption key={entry.version} value={entry.version}>
          {entry.label}
        </NativeSelectOption>
      ))}
    </NativeSelect>
  );
}

export async function targetUrl(target: string, locale: Locale, page: string): Promise<string> {
  return versionPageUrl(target, await manifestOf(target), locale, page);
}
