import type {
  RunComparison,
  RunRecord,
  SettingDifference,
  SettingsComparison,
  TimetreeEstimates,
  YearDate,
} from "@neherlab/app-contracts";
import { runsAuspice, runsCompare, runsGet, runsList, runsResults } from "@neherlab/app-contracts/client";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { NativeSelect, NativeSelectOption } from "@neherlab/design/ui/native-select";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { Link, useNavigate } from "@tanstack/react-router";
import { lazy, Suspense, useCallback, useMemo } from "react";
import ArrowLeftRight from "~icons/lucide/arrow-left-right";

import { useApi, useApiQueries } from "../api/hooks";
import type { ApiCallContext } from "../api/keys";
import { isNotFound } from "../api/retry";
import { LoadingState, PageShell } from "../components/PageShell";
import { Panel } from "../components/Panel";
import { formatLevel, formatSignedDays, yearDateText } from "../format";
import { useLocale } from "../i18n/locale";
import { formatSignificant } from "../results/estimates";
import { fromJsonFloat, nonFiniteLabel } from "../results/numbers";
import { commandSettings } from "../settings/catalog";
import { defaultText, SettingHelpPopover } from "../settings/help";
import { baseName } from "../settings/inputs";
import type { RunTree } from "./CompareTrees";
import { DateIntervals } from "./DateIntervals";
import { MissingRun } from "./MissingRun";
import { ShiftPlot } from "./ShiftPlot";

const SIDE_LABELS = { first: "A", second: "B" } as const;

const RATE_DIGITS = 3;

const CORRELATION_DIGITS = 3;

const CompareTrees = lazy(async () => ({ default: (await import("./CompareTrees")).CompareTrees }));

export function ComparePage({ first, second }: { first: string; second: string }) {
  const [left, right] = useApiQueries(
    [first, second].map((id) => (context: ApiCallContext) => runsGet({ ...context, path: { id } })),
  );

  if (left?.data === undefined || right?.data === undefined) {
    const error = left?.error ?? right?.error ?? null;

    if (error !== null && isNotFound(error)) {
      return (
        <PageShell>
          <MissingRun />
        </PageShell>
      );
    }

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
          {!comparison.same_analysis && (
            <Alert>
              <AlertTitle>Different analyses</AlertTitle>
              <AlertDescription>
                A is {commandSettings(left.command).title}, B is {commandSettings(right.command).title}. The runs have
                no settings or estimates to compare.
              </AlertDescription>
            </Alert>
          )}
          {timetrees && <TimetreeEstimatesComparison left={left} right={right} comparison={comparison} />}
          <SettingsDifferences left={left} right={right} settings={comparison.settings ?? undefined} />
        </>
      )}
      {left.status === "ok" && right.status === "ok" && <TreeComparison left={left} right={right} />}
    </PageShell>
  );
}

function TreeComparison({ left, right }: { left: RunRecord; right: RunRecord }) {
  const runs = [left, right] as const;

  const results = useApiQueries(
    runs.map((run) => (context: ApiCallContext) => runsResults({ ...context, path: { id: run.id } })),
    { staleTime: Infinity },
  );

  const withoutTree = results.flatMap((result, index) =>
    result.data !== undefined && result.data.tree === undefined ? [runs[index]?.title] : [],
  );

  const withTrees = results.every((result) => result.data?.tree !== undefined);

  const documents = useApiQueries(
    runs.map((run) => (context: ApiCallContext) => runsAuspice({ ...context, path: { id: run.id } })),
    { enabled: withTrees, staleTime: Infinity },
  );

  const [firstDocument, secondDocument] = documents.map((query) => query.data);
  const error = [...results, ...documents].find((query) => query.error !== null)?.error ?? null;

  const trees = useMemo<readonly [RunTree, RunTree] | undefined>(
    () =>
      firstDocument === undefined || secondDocument === undefined
        ? undefined
        : [
            { run: figureRun(left), document: firstDocument },
            { run: figureRun(right), document: secondDocument },
          ],
    [firstDocument, left, right, secondDocument],
  );

  if (error !== null) {
    return (
      <Alert variant="destructive">
        <AlertTitle>The trees cannot be read</AlertTitle>
        <AlertDescription>{error.message}</AlertDescription>
      </Alert>
    );
  }

  if (withoutTree.length > 0) {
    return (
      <Alert>
        <AlertTitle>No tree comparison</AlertTitle>
        <AlertDescription>
          {withoutTree.join(" and ")} wrote no Auspice tree, so the trees are not drawn.
        </AlertDescription>
      </Alert>
    );
  }

  if (trees === undefined) {
    return <LoadingState text="Loading the trees" />;
  }

  return (
    <Panel title="Trees" hint="Lines join the tips of the same name in the two trees">
      <Suspense fallback={<LoadingState text="Loading the tree view" />}>
        <CompareTrees first={trees[0]} second={trees[1]} />
      </Suspense>
    </Panel>
  );
}

function figureRun(record: RunRecord) {
  return { id: record.id, title: record.title, command: record.command };
}

function RunPicker({ current, other, side }: { current: RunRecord; other: RunRecord; side: Side }) {
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
    <span className="inline-flex items-center gap-1.5">
      <span className="text-muted-foreground font-bold">{SIDE_LABELS[side]}</span>
      <NativeSelect
        size="sm"
        aria-label={side === "first" ? "Run A" : "Run B"}
        value={current.id}
        onChange={onChange}
        className="max-w-[min(18rem,100%)] font-bold"
      >
        {runs.map((run) => (
          <NativeSelectOption key={run.id} value={run.id}>
            {run.title} ({commandSettings(run.command).title})
          </NativeSelectOption>
        ))}
      </NativeSelect>
    </span>
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
            ? settings.same_config
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
                <RunLink side="first" record={left} />
              </TableHead>
              <TableHead>
                <RunLink side="second" record={right} />
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
  const locale = useLocale();
  const spec = commandSettings(record.command).settings.find((candidate) => candidate.key === difference.key);

  return (
    <TableRow>
      <TableCell className="whitespace-normal">
        {spec?.label ?? difference.key} {spec !== undefined && <SettingHelpPopover spec={spec} value={undefined} />}
        {difference.kind === "input" && (
          <span className="text-muted-foreground block text-xs">Different file contents</span>
        )}
      </TableCell>
      <TableCell className="font-mono text-xs break-all whitespace-normal">
        {difference.kind === "input"
          ? difference.first.map(baseName).join(", ")
          : defaultText(difference.first, locale, spec?.options)}
      </TableCell>
      <TableCell className="font-mono text-xs break-all whitespace-normal">
        {difference.kind === "input"
          ? difference.second.map(baseName).join(", ")
          : defaultText(difference.second, locale, spec?.options)}
      </TableCell>
    </TableRow>
  );
}

function RunLink({ side, record }: { side: Side; record: RunRecord }) {
  return (
    <>
      <span className="text-muted-foreground mr-1.5">{SIDE_LABELS[side]}</span>
      <Link
        to="/runs/$id/results"
        params={{ id: record.id }}
        className="text-primary font-bold underline-offset-4 hover:underline"
      >
        {record.title}
      </Link>
    </>
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
                  <RunLink side="first" record={left} />
                </TableHead>
                <TableHead>
                  <RunLink side="second" record={right} />
                </TableHead>
                <TableHead>Difference (B - A)</TableHead>
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
                first={a.r === undefined ? "-" : formatSignificant(a.r, CORRELATION_DIGITS)}
                second={b.r === undefined ? "-" : formatSignificant(b.r, CORRELATION_DIGITS)}
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
                  estimates.log_likelihood_not_comparable === undefined
                    ? likelihoodChangeText(likelihoodChange)
                    : `not comparable: ${estimates.log_likelihood_not_comparable}`
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
  return date === undefined ? "not dated" : yearDateText(date);
}

function daysText(days: number | undefined): string {
  return days === undefined ? "no interval" : `${Math.round(days)} days`;
}

function daysDifference(days: number | undefined): string {
  return days === undefined ? "-" : formatSignedDays(days);
}

function rateText(estimates: TimetreeEstimates): string {
  const rate = estimates.clock_rate;

  return rate === undefined
    ? "-"
    : `${formatSignificant(rate, RATE_DIGITS, "scientific")}${estimates.clock_rate_fixed ? " (fixed)" : ""}`;
}

function likelihoodText(estimates: TimetreeEstimates): string {
  const written = estimates.log_likelihood;

  if (written === undefined) {
    return "not written";
  }

  const value = fromJsonFloat(written);

  return Number.isFinite(value) ? value.toFixed(1) : `not finite (${nonFiniteLabel(value)})`;
}

function likelihoodChangeText(change: number | undefined): string {
  return change === undefined ? "-" : `${change > 0 ? "+" : ""}${change.toFixed(1)}`;
}

function signedCount(count: number): string {
  return count === 0 ? "0" : `${count > 0 ? "+" : ""}${count}`;
}

type Side = keyof typeof SIDE_LABELS;
