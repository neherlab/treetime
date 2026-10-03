import type { UiTheme } from "@neherlab/app-contracts";
import { create } from "zustand";

interface PreferencesState {
  theme: UiTheme | undefined;
  sidebarWidth: number | undefined;
}

interface PreferencesActions {
  setTheme: (theme: UiTheme) => void;
  setSidebarWidth: (width: number) => void;
}

export const usePreferencesStore = create<PreferencesState & PreferencesActions>()((set) => ({
  theme: undefined,
  sidebarWidth: undefined,
  setTheme: (theme) => {
    set({ theme });
  },
  setSidebarWidth: (sidebarWidth) => {
    set({ sidebarWidth });
  },
}));
