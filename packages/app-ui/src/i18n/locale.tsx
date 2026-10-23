import { createContext, useContext, type ReactNode } from "react";

import { baseLocale, locales, type Locale } from "./generated/runtime.js";

export { m } from "./generated/messages.js";

export { baseLocale, locales, type Locale };

const LocaleContext = createContext<Locale>(baseLocale);

export function LocaleProvider({ locale, children }: { locale: Locale; children: ReactNode }) {
  return <LocaleContext value={locale}>{children}</LocaleContext>;
}

export function useLocale(): Locale {
  return useContext(LocaleContext);
}
