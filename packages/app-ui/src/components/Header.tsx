import { Monitor, Moon, Sun, PanelLeftClose, PanelLeft } from "lucide-react";
import { useTheme } from "next-themes";
import { useCallback } from "react";

import { useVersion } from "../hooks";
import { useAppStore } from "../store/app-store";
import { Button, Tooltip } from "../ui";

const THEME_CYCLE = ["system", "dark", "light"] as const;

type ThemeChoice = (typeof THEME_CYCLE)[number];

const THEME_META: Record<ThemeChoice, { label: string; icon: React.ReactNode }> = {
  system: { label: "System theme", icon: <Monitor size={16} /> },
  dark: { label: "Dark mode", icon: <Moon size={16} /> },
  light: { label: "Light mode", icon: <Sun size={16} /> },
};

function nextTheme(current: string | undefined): ThemeChoice {
  const index = THEME_CYCLE.findIndex((choice) => choice === current);

  return THEME_CYCLE[(index + 1) % THEME_CYCLE.length] ?? "system";
}

function themeMeta(theme: string | undefined): { label: string; icon: React.ReactNode } {
  const choice = THEME_CYCLE.find((candidate) => candidate === theme);

  return choice ? THEME_META[choice] : THEME_META.system;
}

export function Header() {
  const { data: version } = useVersion();
  const { theme, setTheme } = useTheme();
  const sidebarCollapsed = useAppStore((s) => s.sidebarCollapsed);
  const toggleSidebar = useAppStore((s) => s.toggleSidebar);

  const meta = themeMeta(theme);
  const toggleTheme = useCallback(() => setTheme(nextTheme(theme)), [setTheme, theme]);

  return (
    <header className="border-line bg-surface-1 flex h-12 shrink-0 items-center gap-3 border-b px-4">
      <Tooltip.Root>
        <Tooltip.Trigger
          render={
            <Button variant="ghost" size="icon" onClick={toggleSidebar}>
              {sidebarCollapsed ? <PanelLeft size={16} /> : <PanelLeftClose size={16} />}
            </Button>
          }
        />
        <Tooltip.Popup>{sidebarCollapsed ? "Expand sidebar" : "Collapse sidebar"}</Tooltip.Popup>
      </Tooltip.Root>

      <div className="flex items-baseline gap-2">
        <h1 className="text-ink text-base font-semibold tracking-tight">TreeTime</h1>
        {version && <span className="text-2xs text-ink-faint font-mono">v{version.version}</span>}
      </div>

      <div className="flex-1" />

      <Tooltip.Root>
        <Tooltip.Trigger
          render={
            <Button variant="ghost" size="icon" onClick={toggleTheme}>
              {meta.icon}
            </Button>
          }
        />
        <Tooltip.Popup>{meta.label}</Tooltip.Popup>
      </Tooltip.Root>
    </header>
  );
}
