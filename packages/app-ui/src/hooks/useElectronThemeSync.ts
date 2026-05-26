import { useEffect } from "react";
import { useTheme } from "next-themes";

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
