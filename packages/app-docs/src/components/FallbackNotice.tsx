import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/app-ui/ui/alert";
import Languages from "~icons/lucide/languages";

export function FallbackNotice({ locale }: { locale: Locale }) {
  return (
    <LocaleProvider locale={locale}>
      <Alert role="note" className="not-prose mb-6" data-pagefind-ignore>
        <Languages aria-hidden />
        <AlertTitle>{m.docs_fallback_title({}, { locale })}</AlertTitle>
        <AlertDescription>{m.docs_fallback_body({}, { locale })}</AlertDescription>
      </Alert>
    </LocaleProvider>
  );
}
