import type { ReactNode } from "react";
import { I18nextProvider } from "react-i18next";
import { Provider } from "react-redux";
import { ThemeProvider } from "styled-components";

import { AUSPICE_I18N } from "./i18n";
import type { AuspiceStore } from "./store";

const SIDEBAR_THEME = {
  background: "var(--muted)",
  color: "var(--foreground)",
  "font-family": "var(--typeface-sans)",
  sidebarBoxShadow: "rgba(0, 0, 0, 0.15)",
  selectedColor: "var(--primary)",
  unselectedColor: "var(--muted-foreground)",
  alternateBackground: "var(--input)",
};

export function AuspiceProviders({ store, children }: { store: AuspiceStore; children: ReactNode }) {
  return (
    <I18nextProvider i18n={AUSPICE_I18N}>
      <ThemeProvider theme={SIDEBAR_THEME}>
        <Provider store={store}>{children}</Provider>
      </ThemeProvider>
    </I18nextProvider>
  );
}
