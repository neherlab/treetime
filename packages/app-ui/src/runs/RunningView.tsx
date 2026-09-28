import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord } from "@neherlab/app-contracts";
import { runsCancel } from "@neherlab/app-contracts/client";
import { CircleCheck } from "lucide-react";
import { useCallback, useMemo, useState } from "react";

import { useApiMutation } from "../api/hooks";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { formatRate } from "../format";
import { nonFiniteLabel } from "../results/numbers";
import { filterLog, type IterationPoint, type RunProgress } from "../results/progress";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Progress, ProgressLabel, ProgressValue } from "../ui/progress";
import { Spinner } from "../ui/spinner";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { useToastManager } from "../ui/toast";
import { LogLines } from "./LogLines";
import { Plate } from "./Plate";
import { RateTrace } from "./RateTrace";

type LiveFilter = "all" | "warnings";

const LIVE_FILTERS: ReadonlyArray<{ value: LiveFilter; label: string }> = [
  { value: "all", label: "All" },
  { value: "warnings", label: "Warnings" },
];

export function RunningView({ record, progress }: { record: RunRecord; progress: RunProgress }) {
  const [filter, setFilter] = useState<LiveFilter>("all");
  const entries = useMemo(() => filterLog(progress.entries, filter, ""), [filter, progress.entries]);
  const percent = Math.round(Math.min(1, Math.max(0, progress.fraction)) * 100);

  return (
    <div className="grid gap-4 @3xl:grid-cols-[minmax(0,1fr)_26rem]">
      <Panel
        title="Log"
        hint="Follows new lines until you scroll up"
        actions={<OptionToggle label="Log filter" value={filter} onChange={setFilter} options={LIVE_FILTERS} />}
      >
        <div className="p-2.5">
          <LogLines entries={entries} className="h-[32rem]" empty="No log lines yet." />
        </div>
      </Panel>
      <div className="grid content-start gap-4">
        <Panel title="Progress" hint={progress.message === "" ? undefined : progress.message}>
          <div className="grid gap-3 p-3.5">
            <Progress value={percent} aria-label="Run progress">
              <ProgressLabel className="sr-only">Run progress</ProgressLabel>
              <ProgressValue />
            </Progress>
            <ol className="grid gap-1.5">
              {progress.stages.map((stage) => (
                <li key={stage.name} className="flex items-center gap-2">
                  {stage.endSeconds === undefined ? (
                    <Spinner aria-label="Running" className="text-primary size-3.5" />
                  ) : (
                    <CircleCheck aria-label="Done" className="text-success size-3.5" />
                  )}
                  <span className={cn(stage.endSeconds === undefined && "font-medium")}>{stage.name}</span>
                  <span className="text-muted-foreground ml-auto text-xs">
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
              <NumberCell value={point.logLhTotal} format={fixed1} />
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Plate>
  );
}

function NumberCell({ value, format }: { value: number | undefined; format: (value: number) => string }) {
  if (value === undefined) {
    return <TableCell className="text-muted-foreground text-right">-</TableCell>;
  }

  return Number.isFinite(value) ? (
    <TableCell className="text-right">{format(value)}</TableCell>
  ) : (
    <TableCell className="text-destructive text-right font-medium" title="Not a finite number">
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
    <Button type="button" variant="destructive" onClick={onCancel} disabled={busy}>
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
