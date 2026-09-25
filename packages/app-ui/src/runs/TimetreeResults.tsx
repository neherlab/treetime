import type { RunRecordResult, RunResultsResult } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { formatDecimalDate, formatDuration, formatRate } from "../format";
import { fromJsonFloat, nonFiniteLabel } from "../results/numbers";
import type { CoalescentPrior, TimetreeData, TimetreeEstimates } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { SummaryStrip, type SummaryEntry } from "./Panel";
import { Plate } from "./Plate";
import { rttLine, rttPoints } from "./rootToTip";
import { RootToTipPlot } from "./RootToTipPlot";
import { SelectionPanel } from "./SelectionPanel";
import { SkylinePlot } from "./SkylinePlot";
import { initialColorBy, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

const TIMETREE_COLORINGS = ["num_date"];

export function TimetreeResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecordResult;
  results: RunResultsResult;
  data: TimetreeData;
  tree: TreeData | undefined;
}) {
  const estimates = data.estimates;

  if (tree === undefined || estimates === null || estimates === undefined) {
    return <MissingTree />;
  }

  return <TimetreeView record={record} results={results} data={data} estimates={estimates} tree={tree} />;
}

function TimetreeView({
  record,
  results,
  data,
  estimates,
  tree,
}: {
  record: RunRecordResult;
  results: RunResultsResult;
  data: TimetreeData;
  estimates: TimetreeEstimates;
  tree: TreeData;
}) {
  const summary = useMemo(() => timetreeSummary(record, estimates), [estimates, record]);

  const aside = useCallback(
    (link: TreeLink) => <TimetreeAside record={record} tree={tree} data={data} link={link} />,
    [data, record, tree],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      <TreeView data={tree} colorBy={initialColorBy(tree.tree, TIMETREE_COLORINGS)} aside={aside} />
      <OutputFiles record={record} methods={results.methods ?? undefined} citation={results.citation} />
    </div>
  );
}

export function MissingTree() {
  return (
    <p className="border-line bg-surface-1 text-ink-muted rounded-lg border px-4 py-3.5">
      This run wrote no Auspice tree, so the tree view is not available. The output files are listed below.
    </p>
  );
}

function TimetreeAside({
  record,
  tree,
  data,
  link,
}: {
  record: RunRecordResult;
  tree: TreeData;
  data: TimetreeData;
  link: TreeLink;
}) {
  const regression = data.root_to_tip;

  const points = useMemo(
    () => (regression === null || regression === undefined ? [] : rttPoints(regression)),
    [regression],
  );

  const line = useMemo(
    () => (regression === null || regression === undefined ? undefined : rttLine(regression)),
    [regression],
  );

  return (
    <>
      <SelectionPanel runId={record.id} title={record.title} tree={tree.tree} link={link} />
      <Plate
        title="Root-to-tip regression"
        caption="Samples as TreeTime's final clock model saw them; red points are clock-filter outliers, grey points have an inferred date"
      >
        {points.length === 0 ? (
          <p className="m-0 px-2 py-3 text-sm text-[#4b5f5a]">This run wrote no clock regression table.</p>
        ) : (
          <RootToTipPlot
            points={points}
            line={line}
            selected={link.focus?.name}
            inView={link.zoomed ? link.inView : undefined}
            onSelect={link.select}
          />
        )}
      </Plate>
      {data.skyline.length > 0 && (
        <Plate
          title="Effective population size"
          caption="Skyline estimate with its interval, from the coalescent table"
        >
          <SkylinePlot segments={data.skyline} />
        </Plate>
      )}
    </>
  );
}

function timetreeSummary(record: RunRecordResult, estimates: TimetreeEstimates): SummaryEntry[] {
  const relax = estimates.relaxed_clock;
  const r = estimates.r ?? undefined;
  const excluded = estimates.excluded_samples;

  return [
    rootEntry(estimates),
    rateEntry(estimates),
    {
      label: "Temporal signal",
      value: r === undefined ? "not computed" : `r = ${r.toFixed(3)}`,
      detail:
        r === undefined
          ? estimates.clock_rate_fixed
            ? "Fixed clock rate"
            : undefined
          : `R² = ${(estimates.r_squared ?? r * r).toFixed(3)}`,
    },
    {
      label: "Samples",
      value: String(estimates.samples),
      detail: excluded === 0 ? "All in the clock model" : `${excluded} without usable date or clock outliers`,
      tone: excluded === 0 ? undefined : "caution",
    },
    {
      label: "Coalescent prior",
      value: coalescentPriorText(estimates.coalescent_prior),
      detail:
        relax === null || relax === undefined
          ? "Strict clock"
          : `Relaxed clock: slack ${relax.slack}, coupling ${relax.coupling}`,
    },
    likelihoodEntry(record, estimates),
  ];
}

function coalescentPriorText(prior: CoalescentPrior): string {
  if (prior.kind === "fixed") {
    return `Constant size, Tc = ${prior.tc} years`;
  }

  if (prior.kind === "optimized") {
    return "Constant size, optimized Tc";
  }

  if (prior.kind === "skyline") {
    return `Skyline, ${prior.points} points, stiffness ${prior.stiffness}`;
  }

  return "None";
}

function rootEntry(estimates: TimetreeEstimates): SummaryEntry {
  const date = estimates.root_date;

  if (date === null || date === undefined) {
    return { label: "Root date", value: "not dated" };
  }

  const interval = estimates.root_interval;
  const edge = estimates.root_near_interval_edge;

  return {
    label: "Root date",
    value: formatDecimalDate(date),
    detail:
      interval === null || interval === undefined
        ? "No interval computed"
        : `90%: ${formatDecimalDate(interval.lower)} to ${formatDecimalDate(interval.upper)}${edge ? "; the estimate lies at the interval edge" : ""}`,
    tone: edge ? "caution" : undefined,
  };
}

function rateEntry(estimates: TimetreeEstimates): SummaryEntry {
  const rate = estimates.clock_rate;

  if (rate === null || rate === undefined) {
    return { label: "Clock rate", value: "not written" };
  }

  const std = estimates.clock_rate_std ?? undefined;

  const detail = estimates.clock_rate_fixed
    ? std === undefined
      ? "Fixed"
      : `Fixed, std. dev. ${formatRate(std)}`
    : std === undefined
      ? "No standard deviation computed"
      : `± ${formatRate(std)} (1 std. dev.)`;

  return { label: "Clock rate", value: `${formatRate(rate)} /site/yr`, detail };
}

function likelihoodEntry(record: RunRecordResult, estimates: TimetreeEstimates): SummaryEntry {
  const written = estimates.log_likelihood;

  const runTime =
    record.duration_seconds === null || record.duration_seconds === undefined
      ? ""
      : `, ${formatDuration(record.duration_seconds)}`;

  const iterations = `${estimates.iterations} ${estimates.iterations === 1 ? "iteration" : "iterations"}${runTime}`;

  if (written === null || written === undefined) {
    return { label: "Log likelihood", value: "not written", detail: iterations };
  }

  const value = fromJsonFloat(written);

  return Number.isFinite(value)
    ? { label: "Log likelihood", value: value.toFixed(1), detail: iterations }
    : {
        label: "Log likelihood",
        value: "not finite",
        detail: `${nonFiniteLabel(value)}; ${iterations}`,
        tone: "fault",
      };
}
