import type { RunRecordResult } from "@neherlab/app-contracts";
import { useQueryClient } from "@tanstack/react-query";
import { CircleCheck, LoaderCircle } from "lucide-react";
import { useCallback, useMemo, useState } from "react";

import { useBridge } from "../BridgeContext";
import { formatRate } from "../format";
import { RUNS_KEY } from "../queries";
import { nonFiniteLabel } from "../results/numbers";
import { filterLog, type IterationPoint, type RunProgress } from "../results/progress";
import { Button, Segmented, Toast, cn } from "../ui";
import { LogLines } from "./LogLines";
import { Panel } from "./Panel";
import { Plate } from "./Plate";
import { RateTrace } from "./RateTrace";

type LiveFilter = "all" | "warnings";

const LIVE_FILTERS: ReadonlyArray<{ value: LiveFilter; label: string }> = [
  { value: "all", label: "All" },
  { value: "warnings", label: "Warnings" },
];

export function RunningView({ record, progress }: { record: RunRecordResult; progress: RunProgress }) {
  const [filter, setFilter] = useState<LiveFilter>("all");
  const entries = useMemo(() => filterLog(progress.entries, filter, ""), [filter, progress.entries]);
  const percent = Math.round(Math.min(1, Math.max(0, progress.fraction)) * 100);

  return (
    <div className="grid gap-3.5 lg:grid-cols-[minmax(0,1fr)_26rem]">
      <Panel
        title="Log"
        hint="Follows new lines until you scroll up"
        actions={<Segmented label="Log filter" value={filter} onChange={setFilter} options={LIVE_FILTERS} />}
      >
        <div className="p-2.5">
          <LogLines entries={entries} className="h-[32rem]" empty="No log lines yet." />
        </div>
      </Panel>
      <div className="grid content-start gap-3.5">
        <Panel title="Progress" hint={progress.message === "" ? undefined : progress.message}>
          <div className="grid gap-2.5 px-3.5 py-3">
            <progress
              value={percent}
              max={100}
              aria-label="Run progress"
              className="bg-surface-3 [&::-moz-progress-bar]:bg-accent [&::-webkit-progress-bar]:bg-surface-3 [&::-webkit-progress-value]:bg-accent h-2 w-full appearance-none overflow-hidden rounded-sm border-0"
            />
            <ol className="m-0 grid list-none gap-1 p-0">
              {progress.stages.map((stage) => (
                <li key={stage.name} className="flex items-center gap-2">
                  {stage.endSeconds === undefined ? (
                    <LoaderCircle size={13} aria-label="Running" className="text-accent animate-spin" />
                  ) : (
                    <CircleCheck size={13} aria-label="Done" className="text-signal-ok" />
                  )}
                  <span className={cn(stage.endSeconds === undefined && "font-bold")}>{stage.name}</span>
                  <span className="text-ink-faint ml-auto text-xs tabular-nums">
                    {((stage.endSeconds ?? stage.startSeconds) - stage.startSeconds).toFixed(1)} s
                  </span>
                </li>
              ))}
            </ol>
          </div>
        </Panel>
        {progress.iterations.length > 0 && <IterationPanel iterations={progress.iterations} />}
        {record.status === "running" && <CancelButton id={record.id} />}
      </div>
    </div>
  );
}

function IterationPanel({ iterations }: { iterations: readonly IterationPoint[] }) {
  return (
    <Plate title="Clock rate by iteration" caption={latestCaption(iterations.at(-1))}>
      <RateTrace iterations={iterations} />
      <table className="mt-1.5 w-full border-collapse text-right text-xs tabular-nums">
        <thead>
          <tr className="text-[#4b5f5a]">
            <th className="px-1.5 py-1 text-left font-normal">Iteration</th>
            <th className="px-1.5 py-1 font-normal">Rate</th>
            <th className="px-1.5 py-1 font-normal">max Δt</th>
            <th className="px-1.5 py-1 font-normal">rms Δt</th>
            <th className="px-1.5 py-1 font-normal">log L</th>
          </tr>
        </thead>
        <tbody>
          {iterations.map((point) => (
            <tr key={point.iteration} className="border-t border-[#e3e9e7]">
              <td className="px-1.5 py-0.5 text-left">{point.iteration}</td>
              <NumberCell value={point.clockRate} format={formatRate} />
              <NumberCell value={point.maxTimeChange} format={fixed4} />
              <NumberCell value={point.rmsTimeChange} format={fixed4} />
              <NumberCell value={point.logLhTotal} format={fixed1} />
            </tr>
          ))}
        </tbody>
      </table>
    </Plate>
  );
}

function NumberCell({ value, format }: { value: number | undefined; format: (value: number) => string }) {
  if (value === undefined) {
    return <td className="px-1.5 py-0.5 text-[#8a9a95]">-</td>;
  }

  return Number.isFinite(value) ? (
    <td className="px-1.5 py-0.5">{format(value)}</td>
  ) : (
    <td className="px-1.5 py-0.5 font-bold text-[#b42318]" title="Not a finite number">
      {nonFiniteLabel(value)}
    </td>
  );
}

function CancelButton({ id }: { id: string }) {
  const bridge = useBridge();
  const queryClient = useQueryClient();
  const toasts = Toast.useToastManager();
  const [busy, setBusy] = useState(false);

  const cancel = useCallback(async () => {
    setBusy(true);

    try {
      await bridge.cancelRun(id);
      await queryClient.invalidateQueries({ queryKey: RUNS_KEY });
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be cancelled", description: error instanceof Error ? error.message : "" });
    } finally {
      setBusy(false);
    }
  }, [bridge, id, queryClient, toasts]);

  const onCancel = useCallback(() => void cancel(), [cancel]);

  return (
    <Button type="button" variant="outline" onClick={onCancel} disabled={busy}>
      Cancel run
    </Button>
  );
}

function latestCaption(last: IterationPoint | undefined): string | undefined {
  if (last === undefined) {
    return undefined;
  }

  const rate = `Latest ${formatNumber(last.clockRate, formatRate)}`;

  return last.rSquared === undefined ? rate : `${rate}, R² ${formatNumber(last.rSquared, fixed3)}`;
}

function formatNumber(value: number, format: (value: number) => string): string {
  return Number.isFinite(value) ? format(value) : nonFiniteLabel(value);
}

function fixed4(value: number): string {
  return value.toFixed(4);
}

function fixed3(value: number): string {
  return value.toFixed(3);
}

function fixed1(value: number): string {
  return value.toFixed(1);
}
