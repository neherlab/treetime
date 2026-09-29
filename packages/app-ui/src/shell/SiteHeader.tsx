import { version as getVersion } from "@neherlab/app-contracts/client";
import { formatForDisplay } from "@tanstack/react-hotkeys";
import { useTheme } from "next-themes";
import { useCallback, type ComponentType, type SVGProps } from "react";
import Monitor from "~icons/lucide/monitor";
import Moon from "~icons/lucide/moon";
import Search from "~icons/lucide/search";
import Sun from "~icons/lucide/sun";

import { useApi } from "../api/hooks";
import { PALETTE_HOTKEY } from "../hotkeys";
import { useShellStore } from "../store/shell";
import { Button } from "../ui/button";
import { Kbd } from "../ui/kbd";
import { Separator } from "../ui/separator";
import { SidebarTrigger } from "../ui/sidebar";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const THEME_CYCLE = ["system", "light", "dark"] as const;

type ThemeChoice = (typeof THEME_CYCLE)[number];

const THEME_META: Record<ThemeChoice, { label: string; icon: ComponentType<SVGProps<SVGSVGElement>> }> = {
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
      <div className="flex items-center gap-2 px-1">
        <img src={`${import.meta.env.BASE_URL}logo-small.svg`} alt="" className="size-6" />
        <span className="font-heading text-base font-bold">TreeTime</span>
        {version !== undefined && <span className="text-muted-foreground font-mono text-xs">v{version.version}</span>}
      </div>
      <Button variant="outline" onClick={openPalette} className="text-muted-foreground ml-auto justify-start sm:w-72">
        <Search aria-hidden />
        <span className="hidden flex-1 text-left sm:inline">Search runs, settings, examples</span>
        <Kbd className="hidden sm:inline-flex">{formatForDisplay(PALETTE_HOTKEY)}</Kbd>
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
