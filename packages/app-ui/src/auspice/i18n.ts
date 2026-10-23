import sidebar from "auspice/src/locales/en/sidebar.json";
import translation from "auspice/src/locales/en/translation.json";
import i18next from "i18next";

import type { AuspiceTranslate } from "./state";

export const AUSPICE_I18N = i18next.createInstance({
  resources: { en: { sidebar, translation } },
  lng: "en",
  fallbackLng: "en",
  defaultNS: "translation",
  interpolation: { escapeValue: false },
  initImmediate: false,
});

await AUSPICE_I18N.init();

export function translateAuspice(...[key, values]: Parameters<AuspiceTranslate>): string {
  return AUSPICE_I18N.t(key, values);
}
