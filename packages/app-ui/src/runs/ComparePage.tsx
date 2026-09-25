import type { RunRecordResult } from "@neherlab/app-contracts";
import { Link, useNavigate } from "@tanstack/react-router";
import { ArrowLeftRight } from "lucide-react";
import { useCallback, useMemo } from "react";

import { defaultText } from "../analysis/SettingField";
import { formatDecimalDate, formatRate, formatSignedDays } from "../format";
import { useRunList, useRunRecords, useRunResults } from "../queries";
import { ancestorShifts, compareEstimates, timetreeEstimates, type TimetreeEstimates } from "../results/estimates";
import type { RunResults } from "../results/load";
import { nonFiniteLabel } from "../results/numbers";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { settingDifferences, type SettingDifference } from "../settings/differences";
import { baseName } from "../settings/inputs";
import { zJsonObject } from "../settings/json";
import { settingLabel } from "../settings/labels";
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
      {left.command === "timetree" && right.command === "timetree" && left.status === "ok" && right.status === "ok" && (
        <TimetreeComparison left={left} right={right} />
      )}
      <SettingsDifferences left={left} right={right} />
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

function SettingsDifferences({ left, right }: { left: RunRecordResult; right: RunRecordResult }) {
  const differences = useMemo(() => {
    if (left.command !== right.command) {
      return undefined;
    }

    return settingDifferences(
      COMMAND_SETTINGS[left.command].specs,
      { config: zJsonObject.parse(left.config), inputs: left.inputs },
      { config: zJsonObject.parse(right.config), inputs: right.inputs },
    );
  }, [left, right]);

  const sameHash =
    left.config_hash !== null && left.config_hash !== undefined && left.config_hash === right.config_hash;

  return (
    <Panel
      title="Settings that differ"
      hint={
        differences === undefined
          ? "The runs use different analyses"
          : differences.length === 0
            ? sameHash
              ? "None: same settings on the same input contents"
              : "None"
            : `${differences.length} of ${COMMAND_SETTINGS[left.command].specs.length}`
      }
    >
      {differences !== undefined && differences.length > 0 && (
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
              <DifferenceRow key={difference.spec.key} difference={difference} />
            ))}
          </tbody>
        </table>
      )}
    </Panel>
  );
}

function DifferenceRow({ difference }: { difference: SettingDifference }) {
  return (
    <tr className="border-line border-t">
      <td className="px-3.5 py-1.5">
        {settingLabel(difference.spec.key)}{" "}
        <code className="text-ink-faint font-mono text-xs">{difference.spec.flag}</code>
        {difference.kind === "input" && (
          <span className="text-ink-muted block text-xs">
            {difference.sameContent ? "Same file contents" : "Different file contents"}
          </span>
        )}
      </td>
      <td className="px-3.5 py-1.5 font-mono text-xs break-all">
        {difference.kind === "input" ? difference.first.map(baseName).join(", ") : defaultText(difference.first)}
      </td>
      <td className="px-3.5 py-1.5 font-mono text-xs break-all">
        {difference.kind === "input" ? difference.second.map(baseName).join(", ") : defaultText(difference.second)}
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

function TimetreeComparison({ left, right }: { left: RunRecordResult; right: RunRecordResult }) {
  const first = useRunResults(left.id, left.command, true);
  const second = useRunResults(right.id, right.command, true);

  if (first.data === undefined || second.data === undefined) {
    const error = first.error ?? second.error;

    return (
      <p className="text-ink-muted">
        {error === null ? "Reading the outputs..." : `The outputs cannot be read: ${error.message}`}
      </p>
    );
  }

  return <TimetreeEstimatesComparison left={left} right={right} first={first.data} second={second.data} />;
}

function TimetreeEstimatesComparison({
  left,
  right,
  first,
  second,
}: {
  left: RunRecordResult;
  right: RunRecordResult;
  first: RunResults;
  second: RunResults;
}) {
  const firstTree = first.auspice?.tree;
  const secondTree = second.auspice?.tree;

  const estimates = useMemo(() => {
    if (firstTree === undefined || secondTree === undefined) {
      return undefined;
    }

    return [
      timetreeEstimates({
        tree: firstTree,
        clockModel: first.clockModel,
        augurClock: first.augurClock,
        trace: first.trace,
      }),
      timetreeEstimates({
        tree: secondTree,
        clockModel: second.clockModel,
        augurClock: second.augurClock,
        trace: second.trace,
      }),
    ] as const;
  }, [first, firstTree, second, secondTree]);

  const rootRows = useMemo(() => {
    const [a, b] = estimates ?? [];

    return a?.rootDate === undefined || b?.rootDate === undefined
      ? undefined
      : [
          { id: left.id, label: left.title, date: a.rootDate, interval: a.rootInterval, current: false },
          { id: right.id, label: right.title, date: b.rootDate, interval: b.rootInterval, current: true },
        ];
  }, [estimates, left, right]);

  const shifts = useMemo(
    () => (firstTree === undefined || secondTree === undefined ? [] : ancestorShifts(firstTree, secondTree)),
    [firstTree, secondTree],
  );

  if (estimates === undefined || firstTree === undefined) {
    return (
      <p className="text-ink-muted">One of the runs wrote no Auspice tree, so their estimates cannot be compared.</p>
    );
  }

  const [a, b] = estimates;
  const comparison = compareEstimates(a, b);

  const ancestors = firstTree.nodes.length - firstTree.tips.length;

  const meanShift =
    shifts.length === 0 ? undefined : shifts.reduce((sum, shift) => sum + Math.abs(shift.shiftDays), 0) / shifts.length;

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
                first={dateText(a.rootDate)}
                second={dateText(b.rootDate)}
                difference={comparison.rootShiftDays === undefined ? "-" : formatSignedDays(comparison.rootShiftDays)}
              />
              <EstimateRow
                label="90% interval width"
                first={daysText(comparison.intervalWidthDays[0])}
                second={daysText(comparison.intervalWidthDays[1])}
                difference={widthDifference(comparison.intervalWidthDays)}
              />
              <EstimateRow
                label="Clock rate"
                first={rateText(a)}
                second={rateText(b)}
                difference={
                  comparison.ratePercentChange === undefined
                    ? "-"
                    : `${comparison.ratePercentChange > 0 ? "+" : ""}${comparison.ratePercentChange.toFixed(1)}%`
                }
              />
              <EstimateRow
                label="r"
                first={a.r === undefined ? "-" : a.r.toFixed(3)}
                second={b.r === undefined ? "-" : b.r.toFixed(3)}
                difference=""
              />
              <EstimateRow
                label="Excluded samples"
                first={String(a.excludedSamples)}
                second={String(b.excludedSamples)}
                difference={signedCount(b.excludedSamples - a.excludedSamples)}
              />
              <EstimateRow
                label="Log likelihood"
                first={likelihoodText(a)}
                second={likelihoodText(b)}
                difference={likelihoodDifference(a, b)}
              />
            </tbody>
          </table>
        </Panel>
        {rootRows !== undefined && (
          <Plate title="Root date" caption="With the 90% interval of each run">
            <DateIntervals rows={rootRows} />
          </Plate>
        )}
      </div>
      <Plate title="How far each shared ancestor moves" caption={shiftCaption(shifts.length, ancestors, meanShift)}>
        {shifts.length === 0 ? (
          <p className="m-0 px-2 py-3 text-sm text-[#4b5f5a]">
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

function dateText(date: number | undefined): string {
  return date === undefined ? "not dated" : formatDecimalDate(date);
}

function daysText(days: number | undefined): string {
  return days === undefined ? "no interval" : `${Math.round(days)} days`;
}

function widthDifference([first, second]: readonly [number | undefined, number | undefined]): string {
  return first === undefined || second === undefined ? "-" : formatSignedDays(second - first);
}

function rateText(estimates: TimetreeEstimates): string {
  return estimates.rate === undefined ? "-" : `${formatRate(estimates.rate)}${estimates.rateFixed ? " (fixed)" : ""}`;
}

function likelihoodText(estimates: TimetreeEstimates): string {
  const value = estimates.logLikelihood;

  if (value === undefined) {
    return "not written";
  }

  return Number.isFinite(value) ? value.toFixed(1) : `not finite (${nonFiniteLabel(value)})`;
}

function likelihoodDifference(first: TimetreeEstimates, second: TimetreeEstimates): string {
  const a = first.logLikelihood;
  const b = second.logLikelihood;

  return a !== undefined && b !== undefined && Number.isFinite(a) && Number.isFinite(b)
    ? `${b - a > 0 ? "+" : ""}${(b - a).toFixed(1)}`
    : "-";
}

function signedCount(count: number): string {
  return count === 0 ? "0" : `${count > 0 ? "+" : ""}${count}`;
}
