import type { CommandOptions } from "@neherlab/app-contracts";
import { Play, Square, Loader2 } from "lucide-react";
import { useCallback } from "react";

import { useBridge } from "../BridgeContext";
import { useActiveCommand } from "../hooks/useActiveCommand";
import { useAppStore } from "../store/app-store";
import type { FileSlotKind } from "../types";
import { COMMANDS } from "../types";
import { Button } from "../ui";

function buildDataPath(files: Partial<Record<FileSlotKind, { name: string }>>, slot: FileSlotKind): string {
  const file = files[slot];

  if (!file) {
    throw new Error(`Missing required file: ${slot}`);
  }

  return `data/${file.name}`;
}

export function RunButton() {
  const bridge = useBridge();
  const activeCommand = useActiveCommand();
  const files = useAppStore((s) => s.files);
  const runStatus = useAppStore((s) => s.runStatus);
  const progress = useAppStore((s) => s.progress);
  const setRunStatus = useAppStore((s) => s.setRunStatus);
  const setProgress = useAppStore((s) => s.setProgress);
  const setShowResults = useAppStore((s) => s.setShowResults);
  const setAbortController = useAppStore((s) => s.setAbortController);
  const cancelRun = useAppStore((s) => s.cancelRun);

  const commandLabel = COMMANDS.find((c) => c.name === activeCommand)?.label ?? activeCommand;

  const runCommand = useCallback(async () => {
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
      const aln = files.alignment ? { input_fastas: [buildDataPath(files, "alignment")] } : {};
      const dates = files.dates ? { dates: buildDataPath(files, "dates") } : {};
      const states = files.states ? { states: buildDataPath(files, "states") } : {};

      switch (activeCommand) {
        case "timetree":
          await bridge.timetree({ tree, outdir, ...dates, ...aln }, options);
          break;
        case "ancestral":
          await bridge.ancestral({ tree, outdir, ...aln }, options);
          break;
        case "clock":
          await bridge.clock({ tree, outdir, dates: dates.dates ?? "" }, options);
          break;
        case "mugration":
          await bridge.mugration({ tree, outdir, states: states.states ?? "", attribute: "country" }, options);
          break;
        case "optimize":
          await bridge.optimize({ tree, outdir, ...aln }, options);
          break;
        case "prune":
          await bridge.prune({ tree, outdir, ...aln }, options);
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

  const handleRun = useCallback(() => {
    void runCommand();
  }, [runCommand]);

  if (runStatus === "running") {
    const percent = Math.round((progress?.fraction ?? 0) * 100);

    return (
      <div className="space-y-2">
        <div className="flex items-center gap-2">
          <progress
            value={percent}
            max={100}
            className="bg-surface-2 [&::-moz-progress-bar]:bg-accent [&::-webkit-progress-bar]:bg-surface-2 [&::-webkit-progress-value]:bg-accent h-2 w-full appearance-none overflow-hidden rounded-full"
          />
          <span className="text-2xs text-ink-muted shrink-0 font-mono tabular-nums">{percent}%</span>
        </div>
        <div className="flex items-center justify-between">
          <span className="text-2xs text-ink-muted flex items-center gap-1.5">
            <Loader2 size={12} className="animate-spin" />
            {progress?.stage ?? "Starting..."}
          </span>
          <Button variant="outline" size="sm" onClick={cancelRun}>
            <Square size={12} />
            Cancel
          </Button>
        </div>
      </div>
    );
  }

  return (
    <Button variant="solid" size="md" onClick={handleRun} className="w-full">
      <Play size={16} />
      Run {commandLabel}
    </Button>
  );
}
