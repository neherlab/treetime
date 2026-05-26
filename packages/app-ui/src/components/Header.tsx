import { Moon, Sun, PanelLeftClose, PanelLeft } from "lucide-react";
import { useVersion } from "../hooks";
import { useAppStore } from "../store/app-store";

export function Header() {
  const { data: version } = useVersion();
  const darkMode = useAppStore((s) => s.darkMode);
  const toggleDarkMode = useAppStore((s) => s.toggleDarkMode);
  const sidebarCollapsed = useAppStore((s) => s.sidebarCollapsed);
  const toggleSidebar = useAppStore((s) => s.toggleSidebar);

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
        onClick={toggleDarkMode}
        className="rounded p-1.5 text-gray-500 hover:bg-gray-100 hover:text-gray-700 dark:text-gray-400 dark:hover:bg-gray-800 dark:hover:text-gray-200"
        title={darkMode ? "Light mode" : "Dark mode"}
      >
        {darkMode ? <Sun size={18} /> : <Moon size={18} />}
      </button>
    </header>
  );
}
