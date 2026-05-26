import { Clock, Dna, GitBranch, MapPin, Scissors, SlidersHorizontal } from "lucide-react";
import { clsx } from "clsx";
import { useAppStore } from "../store/app-store";
import type { CommandName } from "../types";
import { COMMANDS } from "../types";

const COMMAND_ICONS: Record<CommandName, React.ReactNode> = {
  timetree: <Clock size={20} />,
  ancestral: <Dna size={20} />,
  clock: <GitBranch size={20} />,
  mugration: <MapPin size={20} />,
  optimize: <SlidersHorizontal size={20} />,
  prune: <Scissors size={20} />,
};

export function CommandNav() {
  const activeCommand = useAppStore((s) => s.activeCommand);
  const setActiveCommand = useAppStore((s) => s.setActiveCommand);
  const collapsed = useAppStore((s) => s.sidebarCollapsed);

  return (
    <nav
      className={clsx(
        "flex shrink-0 flex-col gap-1 border-r border-gray-200 bg-gray-50 py-3 dark:border-gray-700 dark:bg-gray-900/50",
        collapsed ? "w-14 px-2" : "w-48 px-3",
      )}
    >
      {COMMANDS.map((cmd) => (
        <button
          key={cmd.name}
          type="button"
          onClick={() => setActiveCommand(cmd.name)}
          title={collapsed ? cmd.label : undefined}
          className={clsx(
            "flex items-center gap-2.5 rounded-md px-2.5 py-2 text-sm font-medium transition-colors",
            activeCommand === cmd.name
              ? "bg-[var(--color-accent-subtle)] text-[var(--color-accent)] dark:bg-[oklch(0.3_0.06_250)] dark:text-[oklch(0.75_0.12_250)]"
              : "text-gray-600 hover:bg-gray-100 hover:text-gray-900 dark:text-gray-400 dark:hover:bg-gray-800 dark:hover:text-gray-200",
          )}
        >
          <span className="shrink-0">{COMMAND_ICONS[cmd.name]}</span>
          {!collapsed && <span className="truncate">{cmd.label}</span>}
        </button>
      ))}
    </nav>
  );
}
