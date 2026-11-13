import { useHotkey } from "@tanstack/react-hotkeys";
import { useNavigate, useParams } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";

import { useNewAnalysis } from "../analysis/useDraftSession";
import { PALETTE_HOTKEY } from "../hotkeys";
import { useShellStore } from "../store/shell";
import { groupedRunOrder, groupRuns } from "./runList";
import { useListedRuns } from "./useListedRuns";

export const RUN_FILTER_ID = "run-filter";

export function useGlobalShortcuts() {
  const navigate = useNavigate();
  const newAnalysis = useNewAnalysis();
  const currentId = useParams({ strict: false, select: (params) => params.id });
  const runFilter = useShellStore((state) => state.runFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);
  const { runs: listed } = useListedRuns(runFilter, commandFilter);
  const runs = useMemo(() => groupedRunOrder(groupRuns(listed)), [listed]);

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
  useHotkey("N", newAnalysis);
  useHotkey("/", () => document.getElementById(RUN_FILTER_ID)?.focus());
  useHotkey("[", () => step(-1));
  useHotkey("]", () => step(1));
}
