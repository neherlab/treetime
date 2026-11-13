import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord } from "@neherlab/app-contracts";
import { runsCancel } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Progress, ProgressLabel, ProgressValue } from "@neherlab/design/ui/progress";
import { Spinner } from "@neherlab/design/ui/spinner";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { useToastManager } from "@neherlab/design/ui/toast";
import { cn } from "cn";
import { useCallback, useMemo, useState } from "react";

import { useApiMutation } from "../api/hooks";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { formatDuration, formatRate } from "../format";
import { nonFiniteLabel } from "../results/numbers";
import { stageSections, type IterationPoint, type LogFilter, type RunProgress } from "../results/progress";
import { RateTrace } from "./RateTrace";
import { STAGE_ICONS, StageLog } from "./StageLog";

const LIVE_FILTERS: ReadonlyArray<{ value: LogFilter; label: string }> = [
  { value: "all", label: "All" },
  { value: "warnings", label: "Warnings" },
];

export function RunningView({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  const [filter, setFilter] = useState<LogFilter>("all");
  const percent = Math.round(Math.min(1, Math.max(0, progress.fraction)) * 100);
  const stages = useMemo(() => stageSections(progress), [progress]);

  return (
    <div className={cn("grid gap-4", progress.iterations.length > 0 && "@3xl:grid-cols-[minmax(0,1fr)_26rem]")}>
      <div className="grid min-w-0 content-start gap-4">
        <Panel
          title="Stages"
          hint={progress.message === "" ? undefined : progress.message}
          actions={record.status === "running" ? <CancelButton id={record.id} /> : undefined}
        >
          <div className="grid gap-3 p-3.5">
            <Progress value={percent} aria-label="Run progress">
              <ProgressLabel className="sr-only">Run progress</ProgressLabel>
              <ProgressValue />
            </Progress>
            {stages.length === 0 ? (
              <p className="text-muted-foreground text-sm">Waiting for the first stage.</p>
            ) : (
              <ol className="grid gap-1.5 text-sm">
                {stages.map((stage) => (
                  <li key={stage.key} className="flex items-center gap-2">
                    {STAGE_ICONS[stage.state]}
                    <span className={cn("min-w-0 flex-1 truncate", stage.state === "running" && "font-bold")}>
                      {stage.name}
                    </span>
                    <span className="text-muted-foreground text-xs tabular-nums">{formatDuration(stage.seconds)}</span>
                  </li>
                ))}
              </ol>
            )}
          </div>
        </Panel>
        <StageLog
          progress={progress}
          filter={filter}
          query=""
          hint="Open stages follow new lines until you scroll up"
          actions={<OptionToggle label="Log filter" value={filter} onChange={setFilter} options={LIVE_FILTERS} />}
        />
      </div>
      {progress.iterations.length > 0 && <IterationPanel iterations={progress.iterations} />}
    </div>
  );
}

function IterationPanel({ iterations }: { iterations: readonly IterationPoint[] }) {
  return (
    <Panel figure title="Clock rate by iteration" hint={latestCaption(iterations.at(-1))}>
      <RateTrace iterations={iterations} />
      <Table className="text-xs">
        <TableHeader>
          <TableRow>
            <TableHead>Iteration</TableHead>
            <TableHead className="text-right">Rate</TableHead>
            <TableHead className="text-right">max Δt</TableHead>
            <TableHead className="text-right">rms Δt</TableHead>
            <TableHead className="text-right">log L</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {iterations.map((point) => (
            <TableRow key={point.iteration}>
              <TableCell>{point.iteration}</TableCell>
              <NumberCell value={point.clockRate} format={formatRate} />
              <NumberCell value={point.maxTimeChange} format={fixed4} />
              <NumberCell value={point.rmsTimeChange} format={fixed4} />
              <NumberCell value={point.logLhTime} format={fixed1} />
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Panel>
  );
}

function NumberCell({ value, format }: { value: number | undefined; format: (value: number) => string }) {
  if (value === undefined) {
    return <TableCell className="text-muted-foreground text-right">-</TableCell>;
  }

  return Number.isFinite(value) ? (
    <TableCell className="text-right">{format(value)}</TableCell>
  ) : (
    <TableCell className="text-destructive text-right font-bold" title="Not a finite number">
      {nonFiniteLabel(value)}
    </TableCell>
  );
}

function CancelButton({ id }: { id: string }) {
  const { mutateAsync: cancelRun } = useApiMutation((context, run: string) =>
    runsCancel({ ...context, path: { id: run } }),
  );

  const toasts = useToastManager();
  const [busy, setBusy] = useState(false);

  const cancel = useCallback(async () => {
    setBusy(true);

    try {
      await cancelRun(id);
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be cancelled", description: errorMessage(error) });
    } finally {
      setBusy(false);
    }
  }, [cancelRun, id, toasts]);

  const onCancel = useCallback(() => void cancel(), [cancel]);

  return (
    <Button type="button" variant="destructive" size="sm" onClick={onCancel} disabled={busy}>
      {busy && <Spinner />}
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
