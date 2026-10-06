import type { MultiplicityRow, SiteHitsRow } from "@neherlab/app-contracts";
import { memo, useMemo } from "react";
import { Bar, BarChart, CartesianGrid, ComposedChart, Label, Line, XAxis, YAxis, type ActiveLabel } from "recharts";

import { ChartContainer, ChartTooltip, ChartTooltipFrame } from "../ui/chart";
import { logAxis, multiplicityPoints, siteHitsPoints } from "./homoplasy";
import { CHART, PLOT_MARGIN, TICK_STYLE } from "./palette";

const EXPECTED_DOT = { r: 2.5, fill: CHART.ink, strokeWidth: 0 };

export const SiteHitsChart = memo(function SiteHitsChart({ rows }: { rows: readonly SiteHitsRow[] }) {
  const points = useMemo(() => siteHitsPoints(rows), [rows]);
  const byHits = useMemo(() => new Map(rows.map((row) => [row.hits, row])), [rows]);
  const yAxis = useMemo(() => logAxis(rows.flatMap((row) => [row.sites, row.expected])), [rows]);

  return (
    <ChartContainer className="aspect-auto h-[260px] w-full">
      <ComposedChart data={points} margin={PLOT_MARGIN} title="Substitutions per site">
        <CartesianGrid vertical={false} stroke={CHART.grid} />
        <XAxis dataKey="hits" tick={TICK_STYLE} stroke={CHART.muted}>
          <Label
            value="Substitutions at a site"
            position="bottom"
            offset={4}
            {...TICK_STYLE}
            className="fill-muted-foreground"
          />
        </XAxis>
        <YAxis
          scale="log"
          domain={yAxis.domain}
          ticks={yAxis.ticks}
          allowDataOverflow
          tickFormatter={logTick}
          tick={TICK_STYLE}
          stroke={CHART.muted}
          width={56}
        >
          <Label value="Sites" angle={-90} position="insideLeft" {...TICK_STYLE} className="fill-muted-foreground" />
        </YAxis>
        <ChartTooltip content={<SiteHitsTooltip byHits={byHits} />} isAnimationActive={false} />
        <Bar dataKey="sitesShown" name="Observed" fill={CHART.accent} isAnimationActive={false} />
        <Line
          dataKey="expectedShown"
          name="Poisson expectation"
          type="linear"
          stroke={CHART.ink}
          strokeWidth={1.5}
          dot={EXPECTED_DOT}
          isAnimationActive={false}
        />
      </ComposedChart>
    </ChartContainer>
  );
});

export const MultiplicityChart = memo(function MultiplicityChart({ rows }: { rows: readonly MultiplicityRow[] }) {
  const points = useMemo(() => multiplicityPoints(rows), [rows]);
  const byBranches = useMemo(() => new Map(rows.map((row) => [row.branches, row])), [rows]);
  const yAxis = useMemo(() => logAxis(rows.map((row) => row.mutations)), [rows]);

  return (
    <ChartContainer className="aspect-auto h-[260px] w-full">
      <BarChart data={points} margin={PLOT_MARGIN} title="Branches per mutation">
        <CartesianGrid vertical={false} stroke={CHART.grid} />
        <XAxis dataKey="branches" tick={TICK_STYLE} stroke={CHART.muted}>
          <Label value="Branches" position="bottom" offset={4} {...TICK_STYLE} className="fill-muted-foreground" />
        </XAxis>
        <YAxis
          scale="log"
          domain={yAxis.domain}
          ticks={yAxis.ticks}
          allowDataOverflow
          tickFormatter={logTick}
          tick={TICK_STYLE}
          stroke={CHART.muted}
          width={56}
        >
          <Label
            value="Distinct substitutions"
            angle={-90}
            position="insideLeft"
            {...TICK_STYLE}
            className="fill-muted-foreground"
          />
        </YAxis>
        <ChartTooltip content={<MultiplicityTooltip byBranches={byBranches} />} isAnimationActive={false} />
        <Bar dataKey="mutationsShown" name="Distinct substitutions" fill={CHART.accent} isAnimationActive={false} />
      </BarChart>
    </ChartContainer>
  );
});

function SiteHitsTooltip({
  byHits,
  active,
  label,
}: {
  byHits: ReadonlyMap<number, SiteHitsRow>;
  active?: boolean;
  label?: ActiveLabel;
}) {
  const row = byHits.get(Number(label));

  if (active !== true || row === undefined) {
    return null;
  }

  return (
    <ChartTooltipFrame>
      <div className="font-bold">
        {row.hits} {row.hits === 1 ? "substitution" : "substitutions"} at a site
      </div>
      <div>Observed sites {row.sites}</div>
      <div>Poisson expectation {row.expected.toPrecision(4)}</div>
    </ChartTooltipFrame>
  );
}

function MultiplicityTooltip({
  byBranches,
  active,
  label,
}: {
  byBranches: ReadonlyMap<number, MultiplicityRow>;
  active?: boolean;
  label?: ActiveLabel;
}) {
  const row = byBranches.get(Number(label));

  if (active !== true || row === undefined) {
    return null;
  }

  return (
    <ChartTooltipFrame>
      <div className="font-bold">
        On {row.branches} {row.branches === 1 ? "branch" : "branches"}
      </div>
      <div>Distinct substitutions {row.mutations}</div>
    </ChartTooltipFrame>
  );
}

function logTick(value: number): string {
  return value < 1 ? String(value) : value.toLocaleString("en-US");
}
