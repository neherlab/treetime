import type { RunRecordResult } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { formatDecimalDate, formatDuration, formatRate } from "../format";
import { timetreeEstimates, type TimetreeEstimates } from "../results/estimates";
import type { AuspiceOutput, RunResults } from "../results/load";
import { coalescentPrior, relaxedClock, timetreeMethods } from "../results/methods";
import { nonFiniteLabel } from "../results/numbers";
import { leastSquares } from "../results/regression";
import { initialColorBy, type ResultTree } from "../results/tree";
import { isNumber, zJsonObject, type JsonObject } from "../settings/json";
import { OutputFiles } from "./OutputFiles";
import { SummaryStrip, type SummaryEntry } from "./Panel";
import { Plate } from "./Plate";
import { RootToTipPlot, type RttLine, type RttPoint } from "./RootToTipPlot";
import { SelectionPanel } from "./SelectionPanel";
import { SkylinePlot } from "./SkylinePlot";
import { TreeView } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

const TIMETREE_COLORINGS = ["num_date"];

export function TimetreeResults({ record, results }: { record: RunRecordResult; results: RunResults }) {
  const config = useMemo(() => zJsonObject.parse(record.config), [record.config]);
  const auspice = results.auspice;

  const estimates = useMemo(
    () =>
      auspice === undefined
        ? undefined
        : timetreeEstimates({
            tree: auspice.tree,
            clockModel: results.clockModel,
            augurClock: results.augurClock,
            trace: results.trace,
          }),
    [auspice, results],
  );

  if (auspice === undefined || estimates === undefined) {
    return <MissingTree />;
  }

  const methods = timetreeMethods(record.treetime_version, config, estimates);

  return (
    <TimetreeView
      record={record}
      config={config}
      estimates={estimates}
      auspice={auspice}
      results={results}
      methods={methods}
    />
  );
}

function TimetreeView({
  record,
  config,
  estimates,
  auspice,
  results,
  methods,
}: {
  record: RunRecordResult;
  config: JsonObject;
  estimates: TimetreeEstimates;
  auspice: AuspiceOutput;
  results: RunResults;
  methods: string;
}) {
  const summary = useMemo(() => timetreeSummary(record, config, estimates), [config, estimates, record]);

  const aside = useCallback(
    (link: TreeLink) => <TimetreeAside record={record} tree={auspice.tree} results={results} link={link} />,
    [auspice.tree, record, results],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      <TreeView auspice={auspice} colorBy={initialColorBy(auspice.tree, TIMETREE_COLORINGS)} aside={aside} />
      <OutputFiles record={record} methods={methods} />
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
  results,
  link,
}: {
  record: RunRecordResult;
  tree: ResultTree;
  results: RunResults;
  link: TreeLink;
}) {
  const points = useMemo(() => rttPoints(tree), [tree]);
  const line = useMemo(() => leastSquaresLine(points), [points]);

  return (
    <>
      <SelectionPanel runId={record.id} title={record.title} tree={tree} link={link} />
      <Plate title="Root-to-tip distance" caption="Every node of the time tree; click a point to select it in the tree">
        <RootToTipPlot
          points={points}
          line={line}
          selected={link.focus?.name}
          inView={link.zoomed ? link.inView : undefined}
          onSelect={link.select}
        />
      </Plate>
      {results.skyline !== undefined && results.skyline.length > 0 && (
        <Plate
          title="Effective population size"
          caption="Skyline estimate with its interval, from the coalescent table"
        >
          <SkylinePlot segments={results.skyline} />
        </Plate>
      )}
    </>
  );
}

function rttPoints(tree: ResultTree): RttPoint[] {
  return tree.nodes.flatMap((node) =>
    node.date === undefined || node.div === undefined
      ? []
      : [
          {
            name: node.name,
            date: node.date,
            div: node.div,
            tip: node.children.length === 0,
            excluded: node.excluded === true && node.children.length === 0,
          },
        ],
  );
}

function leastSquaresLine(points: readonly RttPoint[]): RttLine | undefined {
  const fit = leastSquares(
    points.flatMap((point) => (point.tip && !point.excluded ? [{ x: point.date, y: point.div }] : [])),
  );

  return fit === undefined
    ? undefined
    : {
        ...fit,
        label: `Least squares through the samples in the clock model: slope ${formatRate(fit.slope)} (a visual guide, not TreeTime's rate estimate)`,
      };
}

function timetreeSummary(record: RunRecordResult, config: JsonObject, estimates: TimetreeEstimates): SummaryEntry[] {
  const relax = relaxedClock(config);

  return [
    rootEntry(estimates),
    rateEntry(config, estimates),
    {
      label: "Temporal signal",
      value: estimates.r === undefined ? "not computed" : `r = ${estimates.r.toFixed(3)}`,
      detail:
        estimates.r === undefined
          ? estimates.rateFixed
            ? "Fixed clock rate"
            : undefined
          : `R² = ${(estimates.r ** 2).toFixed(3)}`,
    },
    {
      label: "Samples",
      value: String(estimates.samples),
      detail:
        estimates.excludedSamples === 0
          ? "All in the clock model"
          : `${estimates.excludedSamples} without usable date or clock outliers`,
      tone: estimates.excludedSamples === 0 ? undefined : "caution",
    },
    {
      label: "Coalescent prior",
      value: coalescentPrior(config),
      detail: relax === undefined ? "Strict clock" : `Relaxed clock: slack ${relax[0]}, coupling ${relax[1]}`,
    },
    likelihoodEntry(record, estimates),
  ];
}

function rootEntry(estimates: TimetreeEstimates): SummaryEntry {
  if (estimates.rootDate === undefined) {
    return { label: "Root date", value: "not dated" };
  }

  const interval = estimates.rootInterval;

  return {
    label: "Root date",
    value: formatDecimalDate(estimates.rootDate),
    detail:
      interval === undefined
        ? "No interval computed"
        : `90%: ${formatDecimalDate(interval[0])} to ${formatDecimalDate(interval[1])}${estimates.rootNearIntervalEdge ? "; the estimate lies at the interval edge" : ""}`,
    tone: estimates.rootNearIntervalEdge ? "caution" : undefined,
  };
}

function rateEntry(config: JsonObject, estimates: TimetreeEstimates): SummaryEntry {
  if (estimates.rate === undefined) {
    return { label: "Clock rate", value: "not written" };
  }

  const std = config["clock_std_dev"];
  const fixedDetail = isNumber(std) ? `Fixed, std. dev. ${formatRate(std)}` : "Fixed";

  const estimatedDetail =
    estimates.rateStd === undefined
      ? "No standard deviation computed"
      : `± ${formatRate(estimates.rateStd)} (1 std. dev.)`;

  return {
    label: "Clock rate",
    value: `${formatRate(estimates.rate)} /site/yr`,
    detail: estimates.rateFixed ? fixedDetail : estimatedDetail,
  };
}

function likelihoodEntry(record: RunRecordResult, estimates: TimetreeEstimates): SummaryEntry {
  const value = estimates.logLikelihood;

  const runTime =
    record.duration_seconds === null || record.duration_seconds === undefined
      ? ""
      : `, ${formatDuration(record.duration_seconds)}`;

  const iterations = `${estimates.iterations} ${estimates.iterations === 1 ? "iteration" : "iterations"}${runTime}`;

  if (value === undefined) {
    return { label: "Log likelihood", value: "not written", detail: iterations };
  }

  return Number.isFinite(value)
    ? { label: "Log likelihood", value: value.toFixed(1), detail: iterations }
    : {
        label: "Log likelihood",
        value: "not finite",
        detail: `${nonFiniteLabel(value)}; ${iterations}`,
        tone: "fault",
      };
}
