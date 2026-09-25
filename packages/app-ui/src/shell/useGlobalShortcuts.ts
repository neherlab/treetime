import { useNavigate } from "@tanstack/react-router";
import { useEffect } from "react";

import { useRunList } from "../queries";
import { useShellStore } from "../store/shell";
import { listedRuns } from "./runList";
import { RUN_FILTER_ID } from "./Sidebar";
import { useCurrentRunId } from "./useCurrentRunId";

const TYPING_TAGS = new Set(["INPUT", "TEXTAREA", "SELECT"]);

export function useGlobalShortcuts() {
  const navigate = useNavigate();
  const currentId = useCurrentRunId();
  const { data } = useRunList();
  const runFilter = useShellStore((state) => state.runFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);

  useEffect(() => {
    const runs = listedRuns(data?.runs ?? [], runFilter, commandFilter);

    function onKeyDown(event: KeyboardEvent) {
      const modifier = event.ctrlKey || event.metaKey;

      if (modifier && event.key.toLowerCase() === "k") {
        event.preventDefault();
        setPaletteOpen(true);

        return;
      }

      const target = event.target instanceof HTMLElement ? event.target : null;
      const typing = target !== null && (TYPING_TAGS.has(target.tagName) || target.isContentEditable);

      if (typing || modifier || event.altKey) {
        return;
      }

      if (event.key === "n") {
        event.preventDefault();
        void navigate({ to: "/new" });
      } else if (event.key === "/") {
        event.preventDefault();
        document.getElementById(RUN_FILTER_ID)?.focus();
      } else if (event.key === "[" || event.key === "]") {
        const index = runs.findIndex((run) => run.id === currentId);
        const step = event.key === "]" ? 1 : -1;
        const next = runs[Math.min(runs.length - 1, Math.max(0, index + step))];

        if (next !== undefined) {
          event.preventDefault();
          void navigate({ to: "/runs/$id/results", params: { id: next.id } });
        }
      }
    }

    window.addEventListener("keydown", onKeyDown);

    return () => window.removeEventListener("keydown", onKeyDown);
  }, [commandFilter, currentId, data, navigate, runFilter, setPaletteOpen]);
}
