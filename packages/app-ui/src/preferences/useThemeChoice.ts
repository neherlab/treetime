import type { UiTheme } from "@neherlab/app-contracts";
import { useTheme } from "next-themes";
import { useCallback } from "react";

import { usePreferencesStore } from "../store/preferences";

export function useThemeChoice() {
  const { theme, setTheme } = useTheme();
  const storeTheme = usePreferencesStore((state) => state.setTheme);

  const chooseTheme = useCallback(
    (choice: UiTheme) => {
      setTheme(choice);
      storeTheme(choice);
    },
    [setTheme, storeTheme],
  );

  return { theme, chooseTheme };
}
