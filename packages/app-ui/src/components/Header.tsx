import { Monitor, Moon, Sun, PanelLeftClose, PanelLeft } from "lucide-react";
import { useTheme } from "next-themes";
import { useVersion } from "../hooks";
import { useAppStore } from "../store/app-store";

const THEME_CYCLE = ["system", "dark", "light"] as const;

function ThemeIcon({ theme }: { theme: string | undefined }) {
  switch (theme) {
    case "dark":
      return <Moon size={18} />;
    case "light":
      return <Sun size={18} />;
    default:
      return <Monitor size={18} />;
  }
}

function themeLabel(theme: string | undefined): string {
  switch (theme) {
    case "dark":
      return "Dark mode";
    case "light":
      return "Light mode";
    default:
      return "System theme";
  }
}

export function Header() {
  const { data: version } = useVersion();
  const { theme, setTheme } = useTheme();
  const sidebarCollapsed = useAppStore((s) => s.sidebarCollapsed);
  const toggleSidebar = useAppStore((s) => s.toggleSidebar);

  function cycleTheme() {
    const current = THEME_CYCLE.indexOf(theme as (typeof THEME_CYCLE)[number]);
    const next = THEME_CYCLE[(current + 1) % THEME_CYCLE.length];
    setTheme(next);
  }

  return (
    <header className="flex h-12 shrink-0 items-center gap-3 border-b border-gray-200 bg-white px-4 dark:border-gray-700 dark:bg-gray-900">
      <button
        type="button"
        onClick={toggleSidebar}
        className="rounded p-1.5 text-gray-500 hover:bg-gray-100 hover:text-gray-700 dark:text-gray-400 dark:hover:bg-gray-800 dark:hover:text-gray-200"
        title={sidebarCollapsed ? "Expand sidebar" : "Collapse sidebar"}
      >
        {sidebarCollapsed ? <PanelLeft size={18} /> : <PanelLeftClose size={18} />}
      </button>

      <div className="flex items-baseline gap-2">
        <h1 className="text-base font-semibold text-gray-900 dark:text-gray-100">TreeTime</h1>
        {version && <span className="text-xs text-gray-400 dark:text-gray-500">v{version.version}</span>}
      </div>

      <div className="flex-1" />

      <button
        type="button"
        onClick={cycleTheme}
        className="rounded p-1.5 text-gray-500 hover:bg-gray-100 hover:text-gray-700 dark:text-gray-400 dark:hover:bg-gray-800 dark:hover:text-gray-200"
        title={themeLabel(theme)}
      >
        <ThemeIcon theme={theme} />
      </button>
    </header>
  );
}
