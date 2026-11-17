import type { RunModels, RunRecord } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@neherlab/design/ui/chart";
import { CellText, DataTable, dataColumns, type DataColumn, type RowSelection } from "@neherlab/design/ui/data-table";
import { cn } from "cn";
import { useCallback, useId, useMemo, useState } from "react";
import { Bar, BarChart, CartesianGrid, LabelList, XAxis, YAxis } from "recharts";
import Download from "~icons/lucide/download";

import { Panel } from "../components/Panel";
import { baseName } from "../settings/inputs";
import {
  gammaText,
  indelText,
  modelColumnsInUse,
  modelRows,
  parameterSharing,
  rateText,
  sharedParamsText,
  type ModelColumnsInUse,
  type ModelRow,
  type ParameterSharing,
} from "./modelRows";
import { useSaveRunOutput } from "./OutputFiles";
import { CHART, PLOT_MARGIN, TICK_STYLE } from "./palette";

const LABELED_STATES = 8;

const BAR_RADIUS: [number, number, number, number] = [4, 4, 0, 0];

const FREQUENCY_DOMAIN: [number, string] = [0, "auto"];

const FREQUENCY_MARGIN = { ...PLOT_MARGIN, top: 20 };

const DOCUMENT_ORDER = [{ id: "label", desc: false }];

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

const modelColumn = dataColumns<ModelRow>();

const LABEL_COLUMN = modelColumn.accessor((row) => row.order, {
  id: "label",
  header: "Partition",
  meta: { width: "2fr", minWidth: 180 },
  cell: ({ row }) => <CellText text={row.original.label} />,
});

const ALIGNMENT_COLUMN = modelColumn.accessor((row) => row.alignment, {
  id: "alignment",
  header: "Alignment",
  meta: { width: "1.5fr", minWidth: 120 },
  sortUndefined: "last",
  cell: ({ getValue }) => <CellText text={getValue() ?? ""} />,
});

const SITES_COLUMN = modelColumn.accessor((row) => row.sites, {
  id: "sites",
  header: "Sites",
  meta: { width: 80, numeric: true },
  sortUndefined: "last",
  cell: ({ getValue }) => getValue()?.toLocaleString() ?? "",
});

const MODEL_COLUMN = modelColumn.accessor((row) => row.model, {
  id: "model",
  header: "Model",
  meta: { width: 90 },
});

const MU_COLUMN = modelColumn.accessor((row) => row.mu, {
  id: "mu",
  header: "μ",
  meta: { width: 90, numeric: true },
  cell: ({ getValue }) => rateText(getValue()),
});

const INDEL_COLUMN = modelColumn.accessor((row) => row.indel?.perSite ?? row.indel?.perAlignment, {
  id: "indel",
  header: "Indel rate",
  meta: { width: 100, numeric: true },
  sortUndefined: "last",
  cell: ({ row }) => {
    const { indel } = row.original;

    return indel === undefined ? "" : <IndelRateText text={indelText(indel)} perSite={indel.perSite} />;
  },
});

const GAMMA_COLUMN = modelColumn.accessor((row) => row.gamma?.alpha, {
  id: "gamma",
  header: "α",
  meta: { width: 100, numeric: true, hint: "Gamma shape α of the rate variation across sites, with K rate categories" },
  sortUndefined: "last",
  cell: ({ row }) => (row.original.gamma === undefined ? "" : gammaText(row.original.gamma)),
});

const SHARED_COLUMN = modelColumn.accessor((row) => sharedParamsText(row.shared), {
  id: "shared",
  header: "Shared parameters",
  meta: { width: "2fr", minWidth: 160 },
  cell: ({ getValue }) => <CellText text={getValue()} />,
});

export function SubstitutionModels({ record, models }: { record: RunRecord; models: RunModels }) {
  const save = useSaveRunOutput();
  const file = record.output_files.find((output) => output.kind === "model");

  const onDownload = useCallback(() => {
    if (file !== undefined) {
      void save({ id: record.id, path: file.path, name: baseName(file.path) }, "The model file cannot be downloaded");
    }
  }, [file, record.id, save]);

  return (
    <Panel
      title={models.models.length === 1 ? "Substitution model" : "Substitution models"}
      hint="Select a model to show its equilibrium frequencies and substitution rates"
      actions={
        file !== undefined && (
          <Button type="button" variant="outline" size="sm" onClick={onDownload}>
            <Download aria-hidden />
            Download model (.json)
          </Button>
        )
      }
    >
      <ModelsView models={models} />
    </Panel>
  );
}

export function ModelsView({ models }: { models: RunModels }) {
  const rows = useMemo(() => modelRows(models), [models]);
  const columns = useMemo(() => modelColumns(modelColumnsInUse(rows)), [rows]);

  const [selectedId, setSelectedId] = useState<string | undefined>(undefined);
  const selected = rows.find((row) => row.id === selectedId) ?? rows[0];
  const selection = useMemo<RowSelection>(() => ({ selected: selected?.id, onChange: setSelectedId }), [selected?.id]);

  const sharing = useMemo(
    () => (selected === undefined ? undefined : parameterSharing(rows, selected)),
    [rows, selected],
  );

  return (
    <>
      <DataTable
        label="Substitution models"
        columns={columns}
        rows={rows}
        rowId={rowKey}
        initialSorting={DOCUMENT_ORDER}
        selection={selection}
      />
      {selected !== undefined && <ModelDetail row={selected} sharing={sharing} />}
    </>
  );
}

function modelColumns(inUse: ModelColumnsInUse): Array<DataColumn<ModelRow>> {
  return modelColumn.columns([
    LABEL_COLUMN,
    ...(inUse.alignment ? [ALIGNMENT_COLUMN] : []),
    ...(inUse.sites ? [SITES_COLUMN] : []),
    MODEL_COLUMN,
    MU_COLUMN,
    ...(inUse.indel ? [INDEL_COLUMN] : []),
    ...(inUse.gamma ? [GAMMA_COLUMN] : []),
    ...(inUse.shared ? [SHARED_COLUMN] : []),
  ]);
}

function rowKey(row: ModelRow): string {
  return row.id;
}

function IndelRateText({ text, perSite }: { text: string; perSite: number | undefined }) {
  return <span title={text}>{perSite === undefined ? text : perSite.toPrecision(3)}</span>;
}

function ModelDetail({ row, sharing }: { row: ModelRow; sharing: ParameterSharing | undefined }) {
  const headingId = useId();
  const { states, pi, rates } = row.source;

  return (
    <section aria-labelledby={headingId} className="grid gap-3 border-t p-3">
      <div className="grid gap-1">
        <h3 id={headingId} className="text-sm font-medium">
          {row.label}
        </h3>
        {sharing !== undefined && <SharingNote sharing={sharing} />}
      </div>
      <div className="grid gap-4 @3xl:grid-cols-2">
        <FrequencyChart states={states} pi={pi} />
        <RateHeatmap states={states} rates={rates} />
      </div>
    </section>
  );
}

function SharingNote({ sharing: { params, partitions, all } }: { sharing: ParameterSharing }) {
  return (
    <p className="text-muted-foreground text-xs">
      {all
        ? `Shared parameters: ${params}, shared by all ${partitions.length} partitions`
        : `Shared parameters: ${params}, shared by ${partitions.length} partitions: ${partitions.join(", ")}`}
    </p>
  );
}

function FrequencyChart({ states, pi }: { states: readonly string[]; pi: readonly number[] }) {
  const captionId = useId();

  const data = useMemo(
    () =>
      states.map((state, index) => {
        const value = pi[index] ?? 0;

        return { state, pi: value, label: value.toFixed(3) };
      }),
    [pi, states],
  );

  return (
    <figure aria-labelledby={captionId} className="grid gap-1">
      <figcaption id={captionId} className="text-sm font-medium">
        Equilibrium frequencies
      </figcaption>
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
  const captionId = useId();
  const max = useMemo(() => Math.max(0, ...rates.flat()), [rates]);
  const labeled = states.length <= LABELED_STATES;

  return (
    <figure aria-labelledby={captionId} className="grid content-start gap-1">
      <figcaption id={captionId} className="text-sm font-medium">
        Substitution rates
      </figcaption>
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
