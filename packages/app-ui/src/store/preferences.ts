import type { UiSettings, UiTheme } from "@neherlab/app-contracts";
import { create } from "zustand";

interface PreferencesState {
  preferences: Pick<UiSettings, "theme" | "sidebar_width">;
}

interface PreferencesActions {
  load: (preferences: Pick<UiSettings, "theme" | "sidebar_width">) => void;
  setTheme: (theme: UiTheme) => void;
  setSidebarWidth: (width: number) => void;
}

export const usePreferencesStore = create<PreferencesState & PreferencesActions>()((set) => ({
  preferences: {},
  load: (preferences) => {
    set({ preferences });
  },
  setTheme: (theme) => {
    set(({ preferences }) => ({ preferences: { ...preferences, theme } }));
  },
  setSidebarWidth: (width) => {
    set(({ preferences }) => ({ preferences: { ...preferences, sidebar_width: width } }));
  },
}));
