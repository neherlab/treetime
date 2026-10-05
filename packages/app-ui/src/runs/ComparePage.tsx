import type {
  RunComparison,
  RunRecord,
  SettingDifference,
  SettingsComparison,
  TimetreeEstimates,
  YearDate,
} from "@neherlab/app-contracts";
import { runsCompare, runsGet, runsList } from "@neherlab/app-contracts/client";
import { Link, useNavigate } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";
import ArrowLeftRight from "~icons/lucide/arrow-left-right";

import { defaultText } from "../analysis/SettingField";
import { useApi, useApiQueries } from "../api/hooks";
import type { ApiCallContext } from "../api/keys";
import { LoadingState, PageShell } from "../components/PageShell";
import { Panel } from "../components/Panel";
import { formatLevel, formatRate, formatSignedDays } from "../format";
import { fromJsonFloat, nonFiniteLabel } from "../results/numbers";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { baseName } from "../settings/inputs";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Button } from "../ui/button";
import { NativeSelect, NativeSelectOption } from "../ui/native-select";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { DateIntervals } from "./DateIntervals";
import { ShiftPlot } from "./ShiftPlot";

export function ComparePage({ first, second }: { first: string; second: string }) {
  const [left, right] = useApiQueries(
    [first, second].map((id) => (context: ApiCallContext) => runsGet({ ...context, path: { id } })),
  );

  if (left?.data === undefined || right?.data === undefined) {
    const error = left?.error ?? right?.error ?? null;

    return error === null ? (
      <LoadingState text="Loading the runs" />
    ) : (
      <PageShell>
        <Alert variant="destructive">
          <AlertTitle>The runs cannot be loaded</AlertTitle>
          <AlertDescription>{error.message}</AlertDescription>
        </Alert>
      </PageShell>
    );
  }

  return <Comparison left={left.data} right={right.data} />;
}

function Comparison({ left, right }: { left: RunRecord; right: RunRecord }) {
  const navigate = useNavigate();

  const { data: comparison, error } = useApi(
    (context) => runsCompare({ ...context, path: { id: left.id, other: right.id } }),
    { staleTime: Infinity },
  );

  const timetrees =
    left.command === "timetree" && right.command === "timetree" && left.status === "ok" && right.status === "ok";

  const swap = useCallback(
    () => void navigate({ to: "/compare/$a/$b", params: { a: right.id, b: left.id } }),
    [left.id, navigate, right.id],
  );

  return (
    <PageShell>
      <header className="flex flex-wrap items-center gap-2.5">
        <h1 className="font-heading mr-2 text-2xl font-bold">Compare runs</h1>
        <RunPicker current={left} other={right} side="first" />
        <Button type="button" variant="ghost" size="sm" onClick={swap} aria-label="Swap the runs">
          <ArrowLeftRight aria-hidden />
          Swap
        </Button>
        <RunPicker current={right} other={left} side="second" />
      </header>
      {comparison === undefined &&
        (error === null ? (
          <LoadingState text="Comparing the runs" />
        ) : (
          <Alert variant="destructive">
            <AlertTitle>The runs cannot be compared</AlertTitle>
            <AlertDescription>{error.message}</AlertDescription>
          </Alert>
        ))}
      {comparison !== undefined && (
        <>
          {timetrees && <TimetreeEstimatesComparison left={left} right={right} comparison={comparison} />}
          <SettingsDifferences left={left} right={right} settings={comparison.settings ?? undefined} />
        </>
      )}
    </PageShell>
  );
}

function RunPicker({ current, other, side }: { current: RunRecord; other: RunRecord; side: "first" | "second" }) {
  const { data } = useApi((context) => runsList(context));
  const navigate = useNavigate();
  const runs = useMemo(() => (data?.runs ?? []).filter((run) => run.id !== other.id), [data, other.id]);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => {
      const id = event.target.value;
      const params = side === "first" ? { a: id, b: other.id } : { a: other.id, b: id };

      void navigate({ to: "/compare/$a/$b", params });
    },
    [navigate, other.id, side],
  );

  return (
    <NativeSelect
      size="sm"
      aria-label={side === "first" ? "First run" : "Second run"}
      value={current.id}
      onChange={onChange}
      className="max-w-[min(18rem,100%)] font-bold"
    >
      {runs.map((run) => (
        <NativeSelectOption key={run.id} value={run.id}>
          {run.title} ({COMMAND_INFO[run.command].label})
        </NativeSelectOption>
      ))}
    </NativeSelect>
  );
}

function SettingsDifferences({
  left,
  right,
  settings,
}: {
  left: RunRecord;
  right: RunRecord;
  settings: SettingsComparison | undefined;
}) {
  const differences = settings?.differences ?? [];

  return (
    <Panel
      title="Settings that differ"
      hint={
        settings === undefined
          ? "The runs use different analyses"
          : differences.length === 0
            ? settings.same_config_hash
              ? "None: same settings on the same input contents"
              : "None"
            : `${differences.length} of ${settings.compared}`
      }
    >
      {differences.length > 0 && (
        <Table>
          <TableHeader>
            <TableRow>
              <TableHead>Setting</TableHead>
              <TableHead>
                <RunLink record={left} />
              </TableHead>
              <TableHead>
                <RunLink record={right} />
              </TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {differences.map((difference) => (
              <DifferenceRow key={difference.key} record={left} difference={difference} />
            ))}
          </TableBody>
        </Table>
      )}
    </Panel>
  );
}

function DifferenceRow({ record, difference }: { record: RunRecord; difference: SettingDifference }) {
  const spec = COMMAND_SETTINGS[record.command].settings.find((candidate) => candidate.key === difference.key);

  return (
    <TableRow>
      <TableCell className="whitespace-normal">
        {spec?.label ?? difference.key} <code className="text-muted-foreground font-mono text-xs">{spec?.flag}</code>
        {difference.kind === "input" && (
          <span className="text-muted-foreground block text-xs">
            {difference.same_content ? "Same file contents" : "Different file contents"}
          </span>
        )}
      </TableCell>
      <TableCell className="font-mono text-xs break-all whitespace-normal">
        {difference.kind === "input" ? difference.first.map(baseName).join(", ") : defaultText(difference.first)}
      </TableCell>
      <TableCell className="font-mono text-xs break-all whitespace-normal">
        {difference.kind === "input" ? difference.second.map(baseName).join(", ") : defaultText(difference.second)}
      </TableCell>
    </TableRow>
  );
}

function RunLink({ record }: { record: RunRecord }) {
  return (
    <Link
      to="/runs/$id/results"
      params={{ id: record.id }}
      className="text-primary font-bold underline-offset-4 hover:underline"
    >
      {record.title}
    </Link>
  );
}

function TimetreeEstimatesComparison({
  left,
  right,
  comparison,
}: {
  left: RunRecord;
  right: RunRecord;
  comparison: RunComparison;
}) {
  const estimates = comparison.estimates ?? undefined;
  const ancestors = comparison.ancestors ?? undefined;

  const rootRows = useMemo(() => {
    const a = estimates?.first;
    const b = estimates?.second;

    return a?.root_date === undefined || b?.root_date === undefined
      ? undefined
      : [
          {
            id: left.id,
            label: left.title,
            date: a.root_date,
            interval: a.root_interval ?? undefined,
            current: false,
          },
          {
            id: right.id,
            label: right.title,
            date: b.root_date,
            interval: b.root_interval ?? undefined,
            current: true,
          },
        ];
  }, [estimates, left, right]);

  if (estimates === undefined || ancestors === undefined) {
    return (
      <p className="text-muted-foreground">
        One of the runs wrote no Auspice tree, so their estimates cannot be compared.
      </p>
    );
  }

  const a = estimates.first;
  const b = estimates.second;
  const level = a.root_interval?.level ?? b.root_interval?.level;
  const intervalName = level === undefined ? "Interval" : `${formatLevel(level)} interval`;
  const shifts = ancestors.shifts;
  const rateChange = estimates.clock_rate_change_percent ?? undefined;
  const likelihoodChange = estimates.log_likelihood_change ?? undefined;

  return (
    <>
      <div className="grid gap-4 @4xl:grid-cols-[minmax(0,3fr)_minmax(0,2fr)]">
        <Panel title="Estimates">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead>Estimate</TableHead>
                <TableHead>
                  <RunLink record={left} />
                </TableHead>
                <TableHead>
                  <RunLink record={right} />
                </TableHead>
                <TableHead>Difference</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              <EstimateRow
                label="Root date"
                first={dateText(a.root_date)}
                second={dateText(b.root_date)}
                difference={daysDifference(estimates.root_shift_days)}
              />
              <EstimateRow
                label={`${intervalName} width`}
                first={daysText(a.root_interval?.days)}
                second={daysText(b.root_interval?.days)}
                difference={daysDifference(estimates.root_interval_change_days)}
              />
              <EstimateRow
                label="Clock rate"
                first={rateText(a)}
                second={rateText(b)}
                difference={rateChange === undefined ? "-" : `${rateChange > 0 ? "+" : ""}${rateChange.toFixed(1)}%`}
              />
              <EstimateRow
                label="r"
                first={a.r === undefined ? "-" : a.r.toFixed(3)}
                second={b.r === undefined ? "-" : b.r.toFixed(3)}
                difference=""
              />
              <EstimateRow
                label="Excluded samples"
                first={String(a.excluded_samples)}
                second={String(b.excluded_samples)}
                difference={signedCount(estimates.excluded_samples_change)}
              />
              <EstimateRow
                label="Log likelihood"
                first={likelihoodText(a)}
                second={likelihoodText(b)}
                difference={
                  likelihoodChange === undefined
                    ? "-"
                    : `${likelihoodChange > 0 ? "+" : ""}${likelihoodChange.toFixed(1)}`
                }
              />
            </TableBody>
          </Table>
        </Panel>
        {rootRows !== undefined && (
          <Panel figure title="Root date" hint={`With the ${intervalName.toLowerCase()} of each run`}>
            <DateIntervals rows={rootRows} />
          </Panel>
        )}
      </div>
      <Panel
        figure
        title="How far each shared ancestor moves"
        hint={shiftCaption(shifts.length, ancestors.ancestors, ancestors.mean_absolute_shift_days ?? undefined)}
      >
        {shifts.length === 0 ? (
          <p className="text-muted-foreground px-2 py-3 text-sm">
            The trees share no ancestor with the same set of samples.
          </p>
        ) : (
          <ShiftPlot shifts={shifts} firstLabel={left.title} />
        )}
      </Panel>
    </>
  );
}

function shiftCaption(matched: number, ancestors: number, meanShift: number | undefined): string {
  const counted = `${matched} of ${ancestors} ancestors matched by their set of samples`;

  return meanShift === undefined ? counted : `${counted}; mean absolute shift ${Math.round(meanShift)} days`;
}

function EstimateRow({
  label,
  first,
  second,
  difference,
}: {
  label: string;
  first: string;
  second: string;
  difference: string;
}) {
  return (
    <TableRow>
      <TableCell className="text-muted-foreground">{label}</TableCell>
      <TableCell className="whitespace-normal">{first}</TableCell>
      <TableCell className="whitespace-normal">{second}</TableCell>
      <TableCell className="font-bold">{difference}</TableCell>
    </TableRow>
  );
}

function dateText(date: YearDate | undefined): string {
  return date === undefined ? "not dated" : date.date;
}

function daysText(days: number | undefined): string {
  return days === undefined ? "no interval" : `${Math.round(days)} days`;
}

function daysDifference(days: number | undefined): string {
  return days === undefined ? "-" : formatSignedDays(days);
}

function rateText(estimates: TimetreeEstimates): string {
  const rate = estimates.clock_rate;

  return rate === undefined ? "-" : `${formatRate(rate)}${estimates.clock_rate_fixed ? " (fixed)" : ""}`;
}

function likelihoodText(estimates: TimetreeEstimates): string {
  const written = estimates.log_likelihood;

  if (written === undefined) {
    return "not written";
  }

  const value = fromJsonFloat(written);

  return Number.isFinite(value) ? value.toFixed(1) : `not finite (${nonFiniteLabel(value)})`;
}

function signedCount(count: number): string {
  return count === 0 ? "0" : `${count > 0 ? "+" : ""}${count}`;
}
