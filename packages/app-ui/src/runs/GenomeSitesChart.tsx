import type { HomoplasySite } from "@neherlab/app-contracts";
import { memo, useCallback, useEffect, useMemo, useRef, type RefObject } from "react";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Label,
  Rectangle,
  ReferenceLine,
  XAxis,
  YAxis,
  useActiveTooltipLabel,
  type ActiveLabel,
  type MouseHandlerDataParam,
} from "recharts";

import { ChartContainer, ChartTooltip, ChartTooltipFrame } from "../ui/chart";
import { siteText } from "./homoplasy";
import { CHART, niceAxis, PLOT_MARGIN, THINNED_TICKS, TICK_STYLE } from "./palette";

const MIN_BAR_WIDTH = 1;

const CHART_TITLE = "Sites hit more than once along the genome";

export const GenomeSitesChart = memo(function GenomeSitesChart({
  sites,
  genomeLength,
  zeroBased,
  selected,
  onSelect,
}: {
  sites: readonly HomoplasySite[];
  genomeLength: number;
  zeroBased: boolean;
  selected: number | undefined;
  onSelect: (position: number) => void;
}) {
  const container = useRef<HTMLDivElement>(null);
  const byPosition = useMemo(() => new Map(sites.map((site) => [site.position, site])), [sites]);
  const data = useMemo(() => [...sites], [sites]);
  const domain = useMemo(() => [1, genomeLength], [genomeLength]);

  const ticks = useMemo(
    () => niceAxis(domain).ticks.filter((value) => value >= 1 && value <= genomeLength),
    [domain, genomeLength],
  );

  const offset = zeroBased ? 1 : 0;

  const selectedBar = useMemo(() => {
    const site = selected === undefined ? undefined : byPosition.get(selected);

    return site === undefined ? undefined : barSegment(site);
  }, [byPosition, selected]);

  const tick = useCallback((position: number) => String(position - offset), [offset]);

  const select = useCallback(
    (label: ActiveLabel) => {
      const site = byPosition.get(Number(label));

      if (site !== undefined) {
        onSelect(site.position);
      }
    },
    [byPosition, onSelect],
  );

  const onClick = useCallback((state: MouseHandlerDataParam) => select(state.activeLabel), [select]);

  return (
    <div ref={container}>
      <ChartContainer className="aspect-auto h-[220px] w-full">
        <BarChart data={data} margin={PLOT_MARGIN} onClick={onClick} title={CHART_TITLE} accessibilityLayer>
          <CartesianGrid vertical={false} stroke={CHART.grid} />
          <XAxis
            type="number"
            dataKey="position"
            domain={domain}
            ticks={ticks}
            {...THINNED_TICKS}
            tick={TICK_STYLE}
            tickFormatter={tick}
            stroke={CHART.muted}
          >
            <Label value="Position" position="bottom" offset={4} {...TICK_STYLE} className="fill-muted-foreground" />
          </XAxis>
          <YAxis allowDecimals={false} tick={TICK_STYLE} stroke={CHART.muted} width={48}>
            <Label
              value="Branches"
              angle={-90}
              position="insideLeft"
              {...TICK_STYLE}
              className="fill-muted-foreground"
            />
          </YAxis>
          <ChartTooltip content={<SiteTooltip byPosition={byPosition} />} isAnimationActive={false} />
          <Bar
            dataKey="branches"
            fill={CHART.accent}
            isAnimationActive={false}
            // oxlint-disable-next-line anti-slop/no-shape-in-symbol-names -- `shape` is the Recharts prop that draws each bar
            shape={SiteBar}
          />
          {selectedBar !== undefined && (
            <ReferenceLine segment={selectedBar} stroke={CHART.selection} strokeWidth={MIN_BAR_WIDTH + 1} />
          )}
          <SelectOnEnter container={container} onSelect={select} />
        </BarChart>
      </ChartContainer>
    </div>
  );
});

function SiteBar({ x, y, width, height }: BarGeometry) {
  const drawn = Math.max(MIN_BAR_WIDTH, width);

  return <Rectangle x={x + (width - drawn) / 2} y={y} width={drawn} height={height} fill={CHART.accent} />;
}

interface BarGeometry {
  x: number;
  y: number;
  width: number;
  height: number;
}

function barSegment(site: HomoplasySite) {
  return [
    { x: site.position, y: 0 },
    { x: site.position, y: site.branches },
  ] as const;
}

function SelectOnEnter({
  container,
  onSelect,
}: {
  container: RefObject<HTMLDivElement | null>;
  onSelect: (label: ActiveLabel) => void;
}) {
  const label = useActiveTooltipLabel();

  useEffect(() => {
    const element = container.current;

    const listener = (event: KeyboardEvent) => {
      if (event.key === "Enter") {
        onSelect(label);
      }
    };

    element?.addEventListener("keydown", listener);

    return () => element?.removeEventListener("keydown", listener);
  }, [container, label, onSelect]);

  return null;
}

function SiteTooltip({
  byPosition,
  active,
  label,
}: {
  byPosition: ReadonlyMap<number, HomoplasySite>;
  active?: boolean;
  label?: ActiveLabel;
}) {
  const site = byPosition.get(Number(label));

  if (active !== true || site === undefined) {
    return null;
  }

  return (
    <ChartTooltipFrame>
      <output aria-live="polite" className="font-bold">
        {siteText(site)}
      </output>
      <ul aria-hidden className="grid gap-0.5">
        {site.substitutions.map((substitution) => (
          <li key={substitution.mutation} className="flex justify-between gap-3">
            <span className="font-mono">{substitution.mutation}</span>
            <span className="tabular-nums">{substitution.branches}</span>
          </li>
        ))}
      </ul>
    </ChartTooltipFrame>
  );
}
