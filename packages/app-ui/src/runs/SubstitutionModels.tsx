import type { AlignmentSummary, RunModel, RunModels, RunRecord } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@neherlab/design/ui/chart";
import { cn } from "cn";
import { useCallback, useMemo } from "react";
import { Bar, BarChart, CartesianGrid, LabelList, XAxis, YAxis } from "recharts";
import Download from "~icons/lucide/download";

import { Panel } from "../components/Panel";
import { baseName } from "../settings/inputs";
import { useSaveRunOutput } from "./OutputFiles";
import { CHART, PLOT_MARGIN, TICK_STYLE } from "./palette";

const LABELED_STATES = 8;

const BAR_RADIUS: [number, number, number, number] = [4, 4, 0, 0];

const FREQUENCY_DOMAIN: [number, string] = [0, "auto"];

const FREQUENCY_MARGIN = { ...PLOT_MARGIN, top: 20 };

const RATE_SHADES = [
  "bg-chart-series/10 text-chart-ink",
  "bg-chart-series/20 text-chart-ink",
  "bg-chart-series/30 text-chart-ink",
  "bg-chart-series/40 text-chart-ink",
  "bg-chart-series/50 text-chart-ink",
  "bg-chart-series/60 text-chart-paper",
  "bg-chart-series/70 text-chart-paper",
  "bg-chart-series/80 text-chart-paper",
  "bg-chart-series/90 text-chart-paper",
] as const;

export function SubstitutionModels({
  record,
  models: { alignments, models },
}: {
  record: RunRecord;
  models: RunModels;
}) {
  const save = useSaveRunOutput();
  const file = record.output_files.find((output) => output.kind === "model");

  const onDownload = useCallback(() => {
    if (file !== undefined) {
      void save({ id: record.id, path: file.path, name: baseName(file.path) }, "The model file cannot be downloaded");
    }
  }, [file, record.id, save]);

  return (
    <Panel
      title={models.length === 1 ? "Substitution model" : "Substitution models"}
      hint="Equilibrium frequencies and substitution rates of the fitted model"
      actions={
        file !== undefined && (
          <Button type="button" variant="outline" size="sm" onClick={onDownload}>
            <Download aria-hidden />
            Download model (.json)
          </Button>
        )
      }
    >
      <div className="divide-y">
        {models.map((model) => (
          <ModelView key={modelKey(model)} model={model} alignments={alignments} />
        ))}
      </div>
    </Panel>
  );
}

function ModelView({ model, alignments }: { model: RunModel; alignments: readonly AlignmentSummary[] }) {
  return (
    <section className="grid gap-3 p-3">
      <ModelFacts model={model} alignments={alignments} />
      <div className="grid gap-4 @3xl:grid-cols-2">
        <FrequencyChart states={model.states} pi={model.pi} />
        <RateHeatmap states={model.states} rates={model.rates} />
      </div>
    </section>
  );
}

function ModelFacts({ model, alignments }: { model: RunModel; alignments: readonly AlignmentSummary[] }) {
  const facts =
    model.kind === "sequence"
      ? [
          partitionText(model.partition, model.name),
          model.path === undefined ? `Alignment ${model.alignment}` : baseName(model.path),
          `${model.sites.toLocaleString()} sites`,
          `Model ${model.model}`,
          `μ = ${model.mu.toPrecision(4)}`,
          ...indelText(alignments.find((alignment) => alignment.alignment === model.alignment)),
        ]
      : [`Trait ${model.attribute}`, `Model ${model.model}`, `μ = ${model.mu.toPrecision(4)}`];

  return <p className="text-muted-foreground text-xs">{facts.join(" · ")}</p>;
}

function FrequencyChart({ states, pi }: { states: readonly string[]; pi: readonly number[] }) {
  const data = useMemo(
    () =>
      states.map((state, index) => {
        const value = pi[index] ?? 0;

        return { state, pi: value, label: value.toFixed(3) };
      }),
    [pi, states],
  );

  return (
    <figure className="grid gap-1">
      <figcaption className="text-sm font-medium">Equilibrium frequencies</figcaption>
      <ChartContainer className="aspect-auto h-[200px] w-full">
        <BarChart data={data} margin={FREQUENCY_MARGIN}>
          <CartesianGrid stroke={CHART.grid} vertical={false} />
          <XAxis dataKey="state" tick={TICK_STYLE} interval={0} />
          <YAxis type="number" domain={FREQUENCY_DOMAIN} tick={TICK_STYLE} width={40} />
          <ChartTooltip content={<ChartTooltipContent />} isAnimationActive={false} />
          <Bar dataKey="pi" name="π" fill={CHART.series} radius={BAR_RADIUS} isAnimationActive={false}>
            {states.length <= LABELED_STATES && (
              <LabelList dataKey="label" position="top" fontSize={11} fill={CHART.ink} />
            )}
          </Bar>
        </BarChart>
      </ChartContainer>
    </figure>
  );
}

function RateHeatmap({ states, rates }: { states: readonly string[]; rates: readonly (readonly number[])[] }) {
  const max = useMemo(() => Math.max(0, ...rates.flat()), [rates]);
  const labeled = states.length <= LABELED_STATES;

  return (
    <figure className="grid content-start gap-1">
      <figcaption className="text-sm font-medium">Substitution rates</figcaption>
      <table className="w-full table-fixed border-separate border-spacing-0.5 text-center font-mono text-xs">
        <thead>
          <tr>
            <th className="w-8" aria-label="From state, to state" />
            {states.map((state) => (
              <th key={state} scope="col" className="font-medium">
                {state}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {states.map((from, row) => (
            <tr key={from}>
              <th scope="row" className="font-medium">
                {from}
              </th>
              {states.map((to, column) => (
                <RateCell
                  key={to}
                  from={from}
                  to={to}
                  rate={row === column ? undefined : rates[row]?.[column]}
                  max={max}
                  labeled={labeled}
                />
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <p className="text-muted-foreground text-xs">
        Rate from the row state to the column state per unit of branch length; darker cells change faster.
      </p>
    </figure>
  );
}

function RateCell({
  from,
  to,
  rate,
  max,
  labeled,
}: {
  from: string;
  to: string;
  rate: number | undefined;
  max: number;
  labeled: boolean;
}) {
  if (rate === undefined) {
    return <td className="text-muted-foreground">-</td>;
  }

  const share = max > 0 ? rate / max : 0;
  const shade = RATE_SHADES[Math.min(RATE_SHADES.length - 1, Math.floor(share * RATE_SHADES.length))];

  return (
    <td title={`${from} → ${to}: ${rate.toPrecision(3)}`} className={cn("rounded-sm py-1.5", shade)}>
      {labeled ? rate.toPrecision(2) : <span className="sr-only">{rate.toPrecision(3)}</span>}
    </td>
  );
}

function indelText(alignment: AlignmentSummary | undefined): string[] {
  if (alignment?.indel_rate === undefined) {
    return [];
  }

  const perSite = alignment.sites > 0 ? alignment.indel_rate / alignment.sites : undefined;

  return [
    perSite === undefined
      ? `indel rate ${alignment.indel_rate.toPrecision(4)} per alignment`
      : `indel rate ${perSite.toPrecision(3)} per site (${alignment.indel_rate.toPrecision(4)} per alignment)`,
  ];
}

function partitionText(partition: number, name: string | undefined): string {
  return name === undefined ? `Partition ${partition}` : `Partition ${partition} (${name})`;
}

function modelKey(model: RunModel): string {
  return model.kind === "sequence" ? `sequence-${model.partition}` : `trait-${model.attribute}`;
}
