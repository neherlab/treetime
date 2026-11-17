import type { AlignmentSummary, RunModel, RunModels } from "@neherlab/app-contracts";

import { baseName } from "../settings/inputs";

export function modelRows({ models, alignments }: RunModels): ModelRow[] {
  return models.map((model, order) => ({
    order,
    ...(model.kind === "sequence" ? sequenceRow(model, alignments) : traitRow(model)),
  }));
}

export function modelColumnsInUse(rows: readonly ModelRow[]): ModelColumnsInUse {
  return {
    alignment: rows.some((row) => row.alignment !== undefined),
    sites: rows.some((row) => row.sites !== undefined),
    indel: rows.some((row) => row.indel !== undefined),
    gamma: rows.some((row) => row.gamma !== undefined),
    shared: rows.some((row) => row.shared.length > 0),
  };
}

export function parameterSharing(rows: readonly ModelRow[], selected: ModelRow): ParameterSharing | undefined {
  if (selected.shared.length === 0) {
    return undefined;
  }

  const partitions = rows.flatMap((row) => (row.shared.length > 0 ? [row.label] : []));
  const sequences = rows.filter((row) => row.source.kind === "sequence").length;

  return { params: sharedParamsText(selected.shared), partitions, all: partitions.length === sequences };
}

export function sharedParamsText(shared: readonly SharedParam[]): string {
  return shared.map((param) => SHARED_PARAM_LABELS[param]).join(", ");
}

export function rateText(rate: number): string {
  return rate.toPrecision(4);
}

export function indelText(indel: IndelRate): string {
  return indel.perSite === undefined
    ? `${rateText(indel.perAlignment)} per alignment`
    : `${indel.perSite.toPrecision(3)} per site (${rateText(indel.perAlignment)} per alignment)`;
}

export function gammaText(gamma: RateVariation): string {
  return `${gamma.alpha.toPrecision(3)} (K = ${gamma.categories})`;
}

export interface ModelRow {
  id: string;
  order: number;
  label: string;
  alignment: string | undefined;
  sites: number | undefined;
  model: string;
  mu: number;
  indel: IndelRate | undefined;
  gamma: RateVariation | undefined;
  shared: readonly SharedParam[];
  source: RunModel;
}

export interface IndelRate {
  perSite: number | undefined;
  perAlignment: number;
}

export interface RateVariation {
  alpha: number;
  categories: number;
}

export interface ParameterSharing {
  params: string;
  partitions: string[];
  all: boolean;
}

export interface ModelColumnsInUse {
  alignment: boolean;
  sites: boolean;
  indel: boolean;
  gamma: boolean;
  shared: boolean;
}

export type SharedParam = RunModels["shared_params"][number];

const SHARED_PARAM_LABELS: Record<SharedParam, string> = {
  "exchange-rates": "exchange rates",
  frequencies: "equilibrium frequencies",
  "gamma-shape": "gamma shape",
};

function sequenceRow(model: SequenceModel, alignments: readonly AlignmentSummary[]): Omit<ModelRow, "order"> {
  return {
    id: `sequence-${model.partition}`,
    label: model.name === undefined ? `Partition ${model.partition}` : `Partition ${model.partition} (${model.name})`,
    alignment: model.path === undefined ? `Alignment ${model.alignment}` : baseName(model.path),
    sites: model.sites,
    model: model.model,
    mu: model.mu,
    indel: indelRate(alignments.find((alignment) => alignment.alignment === model.alignment)),
    gamma: model.gamma === undefined ? undefined : { alpha: model.gamma.shape, categories: model.gamma.categories },
    shared: model.shared_params,
    source: model,
  };
}

function traitRow(model: TraitModel): Omit<ModelRow, "order"> {
  return {
    id: `trait-${model.attribute}`,
    label: `Trait ${model.attribute}`,
    alignment: undefined,
    sites: undefined,
    model: model.model,
    mu: model.mu,
    indel: undefined,
    gamma: undefined,
    shared: [],
    source: model,
  };
}

function indelRate(alignment: AlignmentSummary | undefined): IndelRate | undefined {
  if (alignment?.indel_rate === undefined) {
    return undefined;
  }

  return {
    perSite: alignment.sites > 0 ? alignment.indel_rate / alignment.sites : undefined,
    perAlignment: alignment.indel_rate,
  };
}

type SequenceModel = Extract<RunModel, { kind: "sequence" }>;

type TraitModel = Extract<RunModel, { kind: "trait" }>;
