import { runsList } from "@neherlab/app-contracts/client";
import { useHotkey } from "@tanstack/react-hotkeys";
import { useNavigate } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";

import { useApi } from "../api/hooks";
import { PALETTE_HOTKEY } from "../hotkeys";
import { useShellStore } from "../store/shell";
import { listedRuns } from "./runList";
import { useCurrentRunId } from "./useCurrentRunId";

export const RUN_FILTER_ID = "run-filter";

export function useGlobalShortcuts() {
  const navigate = useNavigate();
  const currentId = useCurrentRunId();
  const { data } = useApi((context) => runsList(context));
  const runFilter = useShellStore((state) => state.runFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);
  const runs = useMemo(() => listedRuns(data?.runs ?? [], runFilter, commandFilter), [commandFilter, data, runFilter]);

  const step = useCallback(
    (offset: number) => {
      const index = runs.findIndex((run) => run.id === currentId);
      const next = runs[Math.min(runs.length - 1, Math.max(0, index + offset))];

      if (next !== undefined) {
        void navigate({ to: "/runs/$id/results", params: { id: next.id } });
      }
    },
    [currentId, navigate, runs],
  );

  useHotkey(PALETTE_HOTKEY, () => setPaletteOpen(true));
  useHotkey("N", () => void navigate({ to: "/new" }));
  useHotkey("/", () => document.getElementById(RUN_FILTER_ID)?.focus());
  useHotkey("[", () => step(-1));
  useHotkey("]", () => step(1));
}
