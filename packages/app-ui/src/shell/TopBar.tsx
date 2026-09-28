import { version as getVersion } from "@neherlab/app-contracts/client";
import { Monitor, Moon, Search, Sun } from "lucide-react";
import { useTheme } from "next-themes";
import { useCallback } from "react";

import { useApi } from "../api/hooks";
import { useShellStore } from "../store/shell";
import { Button, Tooltip } from "../ui";

const THEME_CYCLE = ["system", "light", "dark"] as const;

type ThemeChoice = (typeof THEME_CYCLE)[number];

const THEME_META: Record<ThemeChoice, { label: string; icon: React.ReactNode }> = {
  system: { label: "Theme follows the system", icon: <Monitor size={16} /> },
  light: { label: "Light theme", icon: <Sun size={16} /> },
  dark: { label: "Dark theme", icon: <Moon size={16} /> },
};

export function nextTheme(theme: string | undefined): ThemeChoice {
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";

  return THEME_CYCLE[(THEME_CYCLE.indexOf(choice) + 1) % THEME_CYCLE.length] ?? "system";
}

export function TopBar() {
  const { data: version } = useApi((context) => getVersion(context), { staleTime: Infinity });
  const { theme, setTheme } = useTheme();
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";
  const openPalette = useCallback(() => setPaletteOpen(true), [setPaletteOpen]);
  const cycleTheme = useCallback(() => setTheme(nextTheme(theme)), [setTheme, theme]);

  return (
    <header className="border-line bg-surface-1 sticky top-0 z-10 flex items-center gap-3 border-b px-3.5">
      <div className="flex items-center gap-2 text-base font-bold">
        <BrandMark />
        TreeTime
        {version !== undefined && (
          <span className="text-2xs text-ink-faint font-mono font-normal">v{version.version}</span>
        )}
      </div>
      <div className="flex-1" />
      <button
        type="button"
        onClick={openPalette}
        className="border-line bg-surface-2 text-ink-faint hover:border-line-strong flex min-w-0 items-center gap-2.5 rounded-md border py-1 pr-2 pl-2.5 text-left md:min-w-72"
      >
        <Search size={14} aria-hidden />
        <span className="hidden flex-1 md:inline">Search runs, settings, examples</span>
        <kbd className="border-line-strong bg-surface-1 text-ink-muted hidden rounded-sm border px-1 font-mono text-[0.6875rem] md:inline">
          Ctrl K
        </kbd>
      </button>
      <Tooltip.Root>
        <Tooltip.Trigger
          render={
            <Button variant="ghost" size="icon" aria-label={THEME_META[choice].label} onClick={cycleTheme}>
              {THEME_META[choice].icon}
            </Button>
          }
        />
        <Tooltip.Popup>{THEME_META[choice].label}. Click to change.</Tooltip.Popup>
      </Tooltip.Root>
    </header>
  );
}

function BrandMark() {
  return (
    <svg width="22" height="22" viewBox="0 0 22 22" aria-hidden="true" className="text-accent">
      <path
        d="M3 11h4M7 5v12M7 5h5M7 17h7M12 5v-2M12 5v4M12 3h7M12 9h5M14 17v-3M14 17v3M14 14h5M14 20h3"
        fill="none"
        stroke="currentColor"
        strokeWidth="2"
        strokeLinecap="round"
      />
    </svg>
  );
}
