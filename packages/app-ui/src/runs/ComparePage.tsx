import type { RunComparisonResult, RunRecordResult } from "@neherlab/app-contracts";
import { Link, useNavigate } from "@tanstack/react-router";
import { ArrowLeftRight } from "lucide-react";
import { useCallback, useMemo } from "react";

import { defaultText } from "../analysis/SettingField";
import { formatLevel, formatRate, formatSignedDays } from "../format";
import { useRunComparison, useRunList, useRunRecords } from "../queries";
import { fromJsonFloat, nonFiniteLabel } from "../results/numbers";
import type { SettingDifference, SettingsComparison, TimetreeEstimates, YearDate } from "../results/types";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { baseName } from "../settings/inputs";
import { zJsonValue } from "../settings/json";
import { Button } from "../ui";
import { DateIntervals } from "./DateIntervals";
import { Panel } from "./Panel";
import { Plate } from "./Plate";
import { ShiftPlot } from "./ShiftPlot";

export function ComparePage({ first, second }: { first: string; second: string }) {
  const [left, right] = useRunRecords([first, second]);

  if (left?.data === undefined || right?.data === undefined) {
    const error = left?.error ?? right?.error ?? null;

    return (
      <p className="text-ink-muted p-10 text-center">
        {error === null ? "Loading the runs..." : `The runs cannot be loaded: ${error.message}`}
      </p>
    );
  }

  return <Comparison left={left.data} right={right.data} />;
}

function Comparison({ left, right }: { left: RunRecordResult; right: RunRecordResult }) {
  const navigate = useNavigate();
  const { data: comparison, error } = useRunComparison(left.id, right.id, true);

  const timetrees =
    left.command === "timetree" && right.command === "timetree" && left.status === "ok" && right.status === "ok";

  const swap = useCallback(
    () => void navigate({ to: "/compare/$a/$b", params: { a: right.id, b: left.id } }),
    [left.id, navigate, right.id],
  );

  return (
    <div className="mx-auto grid max-w-[110rem] gap-3.5 px-5 pt-4 pb-16">
      <div className="flex flex-wrap items-center gap-2.5">
        <h1 className="mr-2 text-2xl font-bold">Compare runs</h1>
        <RunPicker current={left} other={right} side="first" />
        <Button type="button" variant="ghost" size="sm" onClick={swap} aria-label="Swap the runs">
          <ArrowLeftRight size={14} aria-hidden />
          Swap
        </Button>
        <RunPicker current={right} other={left} side="second" />
      </div>
      {comparison === undefined ? (
        <p className="text-ink-muted">
          {error === null ? "Comparing the runs..." : `The runs cannot be compared: ${error.message}`}
        </p>
      ) : (
        <>
          {timetrees && <TimetreeEstimatesComparison left={left} right={right} comparison={comparison} />}
          <SettingsDifferences left={left} right={right} settings={comparison.settings ?? undefined} />
        </>
      )}
    </div>
  );
}

function RunPicker({
  current,
  other,
  side,
}: {
  current: RunRecordResult;
  other: RunRecordResult;
  side: "first" | "second";
}) {
  const { data } = useRunList();
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
    <select
      aria-label={side === "first" ? "First run" : "Second run"}
      value={current.id}
      onChange={onChange}
      className="border-line-strong bg-surface-1 h-8 max-w-72 rounded-md border px-2 font-bold"
    >
      {runs.map((run) => (
        <option key={run.id} value={run.id}>
          {run.title} ({COMMAND_INFO[run.command].label})
        </option>
      ))}
    </select>
  );
}

function SettingsDifferences({
  left,
  right,
  settings,
}: {
  left: RunRecordResult;
  right: RunRecordResult;
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
        <table className="w-full border-collapse text-left">
          <thead>
            <tr className="text-ink-faint text-xs">
              <th className="px-3.5 py-1.5 font-normal">Setting</th>
              <th className="px-3.5 py-1.5 font-normal">
                <RunLink record={left} />
              </th>
              <th className="px-3.5 py-1.5 font-normal">
                <RunLink record={right} />
              </th>
            </tr>
          </thead>
          <tbody>
            {differences.map((difference) => (
              <DifferenceRow key={difference.key} record={left} difference={difference} />
            ))}
          </tbody>
        </table>
      )}
    </Panel>
  );
}

function DifferenceRow({ record, difference }: { record: RunRecordResult; difference: SettingDifference }) {
  const spec = COMMAND_SETTINGS[record.command].specs.find((candidate) => candidate.key === difference.key);

  return (
    <tr className="border-line border-t">
      <td className="px-3.5 py-1.5">
        {spec?.label ?? difference.key} <code className="text-ink-faint font-mono text-xs">{spec?.flag}</code>
        {difference.kind === "input" && (
          <span className="text-ink-muted block text-xs">
            {difference.same_content ? "Same file contents" : "Different file contents"}
          </span>
        )}
      </td>
      <td className="px-3.5 py-1.5 font-mono text-xs break-all">
        {difference.kind === "input"
          ? difference.first.map(baseName).join(", ")
          : defaultText(zJsonValue.parse(difference.first))}
      </td>
      <td className="px-3.5 py-1.5 font-mono text-xs break-all">
        {difference.kind === "input"
          ? difference.second.map(baseName).join(", ")
          : defaultText(zJsonValue.parse(difference.second))}
      </td>
    </tr>
  );
}

function RunLink({ record }: { record: RunRecordResult }) {
  return (
    <Link to="/runs/$id/results" params={{ id: record.id }} className="text-accent font-bold">
      {record.title}
    </Link>
  );
}

function TimetreeEstimatesComparison({
  left,
  right,
  comparison,
}: {
  left: RunRecordResult;
  right: RunRecordResult;
  comparison: RunComparisonResult;
}) {
  const estimates = comparison.estimates ?? undefined;
  const ancestors = comparison.ancestors ?? undefined;

  const rootRows = useMemo(() => {
    const a = estimates?.first;
    const b = estimates?.second;

    return a?.root_date === null || a?.root_date === undefined || b?.root_date === null || b?.root_date === undefined
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
      <p className="text-ink-muted">One of the runs wrote no Auspice tree, so their estimates cannot be compared.</p>
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
      <div className="grid gap-3.5 xl:grid-cols-[minmax(0,3fr)_minmax(0,2fr)]">
        <Panel title="Estimates">
          <table className="w-full border-collapse text-left text-sm tabular-nums">
            <thead>
              <tr className="text-ink-faint text-xs">
                <th className="px-3.5 py-1.5 font-normal">Estimate</th>
                <th className="px-3.5 py-1.5 font-normal">
                  <RunLink record={left} />
                </th>
                <th className="px-3.5 py-1.5 font-normal">
                  <RunLink record={right} />
                </th>
                <th className="px-3.5 py-1.5 font-normal">Difference</th>
              </tr>
            </thead>
            <tbody>
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
                first={a.r === null || a.r === undefined ? "-" : a.r.toFixed(3)}
                second={b.r === null || b.r === undefined ? "-" : b.r.toFixed(3)}
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
            </tbody>
          </table>
        </Panel>
        {rootRows !== undefined && (
          <Plate title="Root date" caption={`With the ${intervalName.toLowerCase()} of each run`}>
            <DateIntervals rows={rootRows} />
          </Plate>
        )}
      </div>
      <Plate
        title="How far each shared ancestor moves"
        caption={shiftCaption(shifts.length, ancestors.ancestors, ancestors.mean_absolute_shift_days ?? undefined)}
      >
        {shifts.length === 0 ? (
          <p className="text-plate-muted m-0 px-2 py-3 text-sm">
            The trees share no ancestor with the same set of samples.
          </p>
        ) : (
          <ShiftPlot shifts={shifts} firstLabel={left.title} />
        )}
      </Plate>
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
    <tr className="border-line border-t">
      <td className="text-ink-muted px-3.5 py-1.5">{label}</td>
      <td className="px-3.5 py-1.5">{first}</td>
      <td className="px-3.5 py-1.5">{second}</td>
      <td className="px-3.5 py-1.5 font-bold">{difference}</td>
    </tr>
  );
}

function dateText(date: YearDate | null | undefined): string {
  return date === null || date === undefined ? "not dated" : date.date;
}

function daysText(days: number | undefined): string {
  return days === undefined ? "no interval" : `${Math.round(days)} days`;
}

function daysDifference(days: number | null | undefined): string {
  return days === null || days === undefined ? "-" : formatSignedDays(days);
}

function rateText(estimates: TimetreeEstimates): string {
  const rate = estimates.clock_rate;

  return rate === null || rate === undefined
    ? "-"
    : `${formatRate(rate)}${estimates.clock_rate_fixed ? " (fixed)" : ""}`;
}

function likelihoodText(estimates: TimetreeEstimates): string {
  const written = estimates.log_likelihood;

  if (written === null || written === undefined) {
    return "not written";
  }

  const value = fromJsonFloat(written);

  return Number.isFinite(value) ? value.toFixed(1) : `not finite (${nonFiniteLabel(value)})`;
}

function signedCount(count: number): string {
  return count === 0 ? "0" : `${count > 0 ? "+" : ""}${count}`;
}
