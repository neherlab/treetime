import { version as getVersion } from "@neherlab/app-contracts/client";
import { Monitor, Moon, Search, Sun, type LucideIcon } from "lucide-react";
import { useTheme } from "next-themes";
import { useCallback } from "react";

import { useApi } from "../api/hooks";
import { useShellStore } from "../store/shell";
import { Button } from "../ui/button";
import { Kbd, KbdGroup } from "../ui/kbd";
import { Separator } from "../ui/separator";
import { SidebarTrigger } from "../ui/sidebar";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const THEME_CYCLE = ["system", "light", "dark"] as const;

type ThemeChoice = (typeof THEME_CYCLE)[number];

const THEME_META: Record<ThemeChoice, { label: string; icon: LucideIcon }> = {
  system: { label: "Theme follows the system", icon: Monitor },
  light: { label: "Light theme", icon: Sun },
  dark: { label: "Dark theme", icon: Moon },
};

export function nextTheme(theme: string | undefined): ThemeChoice {
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";

  return THEME_CYCLE[(THEME_CYCLE.indexOf(choice) + 1) % THEME_CYCLE.length] ?? "system";
}

export function SiteHeader() {
  const { data: version } = useApi((context) => getVersion(context), { staleTime: Infinity });
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);
  const openPalette = useCallback(() => setPaletteOpen(true), [setPaletteOpen]);

  return (
    <header className="bg-background z-20 flex h-(--header-height) shrink-0 items-center gap-2 border-b px-3">
      <SidebarTrigger />
      <Separator orientation="vertical" className="data-vertical:h-4 data-vertical:self-auto" />
      <div className="flex items-baseline gap-2 px-1">
        <span className="font-heading text-base font-semibold">TreeTime</span>
        {version !== undefined && <span className="text-muted-foreground font-mono text-xs">v{version.version}</span>}
      </div>
      <Button variant="outline" onClick={openPalette} className="text-muted-foreground ml-auto justify-start sm:w-72">
        <Search aria-hidden />
        <span className="hidden flex-1 text-left sm:inline">Search runs, settings, examples</span>
        <KbdGroup className="hidden sm:inline-flex">
          <Kbd>Ctrl</Kbd>
          <Kbd>K</Kbd>
        </KbdGroup>
      </Button>
      <ThemeButton />
    </header>
  );
}

function ThemeButton() {
  const { theme, setTheme } = useTheme();
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";
  const { label, icon: Icon } = THEME_META[choice];
  const cycleTheme = useCallback(() => setTheme(nextTheme(theme)), [setTheme, theme]);

  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button variant="ghost" size="icon" aria-label={label} onClick={cycleTheme}>
            <Icon aria-hidden />
          </Button>
        }
      />
      <TooltipContent>{label}. Click to change.</TooltipContent>
    </Tooltip>
  );
}
