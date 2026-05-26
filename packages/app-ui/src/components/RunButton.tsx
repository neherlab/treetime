import { Play, Square, Loader2 } from "lucide-react";
import { useCallback } from "react";
import { useAppStore } from "../store/app-store";
import { useBridge } from "../BridgeContext";
import type { CommandOptions } from "@neherlab/app-contracts";
import type { FileSlotKind } from "../types";
import { COMMANDS } from "../types";

function buildDataPath(files: Partial<Record<FileSlotKind, { name: string }>>, slot: FileSlotKind): string {
  const file = files[slot];
  if (!file) {
    throw new Error(`Missing required file: ${slot}`);
  }
  return `data/${file.name}`;
}

export function RunButton() {
  const bridge = useBridge();
  const activeCommand = useAppStore((s) => s.activeCommand);
  const files = useAppStore((s) => s.files);
  const runStatus = useAppStore((s) => s.runStatus);
  const progress = useAppStore((s) => s.progress);
  const setRunStatus = useAppStore((s) => s.setRunStatus);
  const setProgress = useAppStore((s) => s.setProgress);
  const setShowResults = useAppStore((s) => s.setShowResults);
  const setAbortController = useAppStore((s) => s.setAbortController);
  const cancelRun = useAppStore((s) => s.cancelRun);

  const commandLabel = COMMANDS.find((c) => c.name === activeCommand)?.label ?? activeCommand;

  const handleRun = useCallback(async () => {
    const controller = new AbortController();
    setAbortController(controller);
    setRunStatus("running");
    setProgress(undefined);

    const options: CommandOptions = {
      signal: controller.signal,
      onProgress(event) {
        setProgress({ stage: event.stage, fraction: event.fraction, message: event.message });
      },
    };

    try {
      const tree = buildDataPath(files, "tree");
      const outdir = `tmp/web/${activeCommand}`;
      const alignment = files.alignment ? buildDataPath(files, "alignment") : undefined;
      const dates = files.dates ? buildDataPath(files, "dates") : undefined;
      const states = files.states ? buildDataPath(files, "states") : undefined;

      switch (activeCommand) {
        case "timetree":
          await bridge.timetree({ tree, dates, input_fastas: alignment ? [alignment] : undefined, outdir }, options);
          break;
        case "ancestral":
          await bridge.ancestral({ tree, input_fastas: alignment ? [alignment] : undefined, outdir }, options);
          break;
        case "clock":
          await bridge.clock({ tree, dates: dates ?? "", outdir }, options);
          break;
        case "mugration":
          await bridge.mugration({ tree, states: states ?? "", attribute: "country", outdir }, options);
          break;
        case "optimize":
          await bridge.optimize({ tree, input_fastas: alignment ? [alignment] : undefined, outdir }, options);
          break;
        case "prune":
          await bridge.prune({ tree, input_fastas: alignment ? [alignment] : undefined, outdir }, options);
          break;
      }
      setRunStatus("completed");
      setShowResults(true);
    } catch {
      if (!controller.signal.aborted) {
        setRunStatus("failed");
      }
    } finally {
      setAbortController(undefined);
    }
  }, [bridge, activeCommand, files, setRunStatus, setProgress, setShowResults, setAbortController]);

  if (runStatus === "running") {
    return (
      <div className="space-y-2">
        <div className="flex items-center gap-2">
          <div className="h-2 flex-1 overflow-hidden rounded-full bg-gray-200 dark:bg-gray-700">
            <div
              className="h-full rounded-full bg-[var(--color-accent)] transition-all duration-300"
              style={{ width: `${(progress?.fraction ?? 0) * 100}%` }}
            />
          </div>
          <span className="shrink-0 text-xs tabular-nums text-gray-500">
            {Math.round((progress?.fraction ?? 0) * 100)}%
          </span>
        </div>
        <div className="flex items-center justify-between">
          <span className="flex items-center gap-1.5 text-xs text-gray-500 dark:text-gray-400">
            <Loader2 size={12} className="animate-spin" />
            {progress?.stage ?? "Starting..."}
          </span>
          <button
            type="button"
            onClick={cancelRun}
            className="flex items-center gap-1 rounded-md border border-gray-300 px-2.5 py-1 text-xs font-medium text-gray-600 hover:bg-gray-50 dark:border-gray-600 dark:text-gray-400 dark:hover:bg-gray-800"
          >
            <Square size={12} />
            Cancel
          </button>
        </div>
      </div>
    );
  }

  return (
    <button
      type="button"
      onClick={handleRun}
      className="flex w-full items-center justify-center gap-2 rounded-lg bg-[var(--color-accent)] px-4 py-2.5 text-sm font-semibold text-white shadow-sm transition-colors hover:bg-[var(--color-accent-hover)]"
    >
      <Play size={16} />
      Run {commandLabel}
    </button>
  );
}
