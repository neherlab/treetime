import { Link } from "@tanstack/react-router";
import { Clock, Dna, GitBranch, MapPin, Scissors, SlidersHorizontal } from "lucide-react";
import { useCallback } from "react";

import { useActiveCommand } from "../hooks/useActiveCommand";
import { useAppStore } from "../store/app-store";
import type { CommandName } from "../types";
import { COMMANDS } from "../types";
import { cn } from "../ui";

const COMMAND_ICONS: Record<CommandName, React.ReactNode> = {
  timetree: <Clock size={18} />,
  ancestral: <Dna size={18} />,
  clock: <GitBranch size={18} />,
  mugration: <MapPin size={18} />,
  optimize: <SlidersHorizontal size={18} />,
  prune: <Scissors size={18} />,
};

export function CommandNav() {
  const activeCommand = useActiveCommand();
  const collapsed = useAppStore((s) => s.sidebarCollapsed);
  const resetRun = useAppStore((s) => s.resetRun);
  const handleNavigate = useCallback(() => resetRun(), [resetRun]);

  return (
    <nav
      className={cn(
        "border-line bg-surface-1 flex shrink-0 flex-col gap-0.5 overflow-y-auto border-r py-3",
        collapsed ? "w-14 px-2" : "w-48 px-3",
      )}
    >
      {COMMANDS.map((cmd) => {
        const active = activeCommand === cmd.name;

        return (
          <Link
            key={cmd.name}
            to="/$command"
            params={{ command: cmd.name }}
            onClick={handleNavigate}
            title={collapsed ? cmd.label : undefined}
            className={cn(
              "flex items-center gap-2.5 rounded-md px-2.5 py-2 text-sm font-medium transition-colors",
              active ? "bg-accent-subtle text-accent" : "text-ink-muted hover:bg-surface-2 hover:text-ink",
            )}
          >
            <span className="shrink-0">{COMMAND_ICONS[cmd.name]}</span>
            {!collapsed && <span className="truncate">{cmd.label}</span>}
          </Link>
        );
      })}
    </nav>
  );
}
