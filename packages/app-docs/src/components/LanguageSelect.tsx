import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { buttonVariants } from "@neherlab/app-ui/ui/button";
import { cn } from "cn";

const LOCALE_NAMES: Record<Locale, (locale: Locale) => string> = {
  en: (locale) => m.docs_locale_en({}, { locale }),
  de: (locale) => m.docs_locale_de({}, { locale }),
};

export function LanguageSelect({ urls, locale }: { urls: ReadonlyMap<Locale, string>; locale: Locale }) {
  return (
    <LocaleProvider locale={locale}>
      <nav aria-label={m.docs_language({}, { locale })} className="flex items-center">
        {[...urls].map(([target, url]) => (
          <a
            key={target}
            href={url}
            hrefLang={target}
            lang={target}
            aria-current={target === locale ? "page" : undefined}
            className={cn(buttonVariants({ variant: "ghost", size: "sm" }), "aria-[current=page]:bg-muted")}
          >
            <span aria-hidden className="uppercase">
              {target}
            </span>
            <span className="sr-only">{LOCALE_NAMES[target](target)}</span>
          </a>
        ))}
      </nav>
    </LocaleProvider>
  );
}
