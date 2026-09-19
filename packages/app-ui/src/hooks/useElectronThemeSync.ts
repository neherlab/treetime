import { useTheme } from "next-themes";
import { useEffect } from "react";

declare global {
  interface Window {
    electronTheme?: { setTheme: (theme: string) => void };
  }
}

export function useElectronThemeSync() {
  const { theme } = useTheme();

  useEffect(() => {
    if (theme) {
      window.electronTheme?.setTheme(theme);
    }
  }, [theme]);
}
