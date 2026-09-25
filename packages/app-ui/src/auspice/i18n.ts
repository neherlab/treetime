import sidebar from "auspice/src/locales/en/sidebar.json";
import translation from "auspice/src/locales/en/translation.json";
import i18next from "i18next";

export const AUSPICE_I18N = i18next.createInstance({
  resources: { en: { sidebar, translation } },
  lng: "en",
  fallbackLng: "en",
  defaultNS: "translation",
  interpolation: { escapeValue: false },
  initImmediate: false,
});

await AUSPICE_I18N.init();
