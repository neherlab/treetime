import { csvParse, tsvParse, type DSVRowString } from "d3-dsv";
import * as z from "zod";

import { parseNumber, parseOptionalNumber } from "./numbers";

export interface ClockModel {
  rate: number;
  intercept: number;
  fixed: boolean;
  r: number | undefined;
}

export interface AugurClock {
  rate: number;
  rateStd: number | undefined;
}

export interface GtrSummary {
  model: string;
  mu: number;
}

export interface TraceRow {
  iteration: number;
  maxTimeChange: number | undefined;
  rmsTimeChange: number | undefined;
  logLhSeq: number | undefined;
  logLhPos: number | undefined;
  logLhCoal: number | undefined;
  logLhTotal: number | undefined;
}

interface Interval {
  value: number;
  lower: number;
  upper: number;
}

export interface SkylineSegment {
  start: number;
  end: number;
  tc: Interval;
  ne: Interval;
}

export interface ClockRow {
  name: string;
  div: number;
  date: number | undefined;
  predictedDate: number;
  deviation: number | undefined;
  outlier: boolean;
}

export interface Traits {
  attribute: string;
  states: ReadonlyMap<string, string>;
}

export function readClockModel(text: string): ClockModel {
  const model = zClockModel.parse(JSON.parse(text));

  return {
    rate: model.clock_rate,
    intercept: model.intercept,
    fixed: model.stats === "fixed",
    r: model.stats === "fixed" ? undefined : model.stats.estimated.r_val,
  };
}

export function readAugurClock(text: string): AugurClock | undefined {
  const clock = zAugurNodeData.parse(JSON.parse(text)).clock;

  return clock === undefined ? undefined : { rate: clock.rate, rateStd: clock.rate_std };
}

export function readTotalBranchLength(text: string): number | undefined {
  const lengths = Object.values(zAugurNodeData.parse(JSON.parse(text)).nodes ?? {}).flatMap(
    (node) => node.branch_length ?? [],
  );

  return lengths.length === 0 ? undefined : lengths.reduce((sum, length) => sum + length, 0);
}

export function readGtr(text: string): GtrSummary {
  const gtr = zGtr.parse(JSON.parse(text));

  return { model: gtr.model_name, mu: gtr.mu };
}

export function readTracelog(text: string): TraceRow[] {
  return csvParse(text).map((row, iteration) => ({
    iteration,
    maxTimeChange: optionalColumn(row, "max_time_change"),
    rmsTimeChange: optionalColumn(row, "rms_time_change"),
    logLhSeq: optionalColumn(row, "log_lh_seq"),
    logLhPos: optionalColumn(row, "log_lh_pos"),
    logLhCoal: optionalColumn(row, "log_lh_coal"),
    logLhTotal: optionalColumn(row, "log_lh_total"),
  }));
}

export function readCoalescentTsv(text: string): SkylineSegment[] {
  return tsvParse(text).map((row) => ({
    start: column(row, "segment.start"),
    end: column(row, "segment.end"),
    tc: { value: column(row, "T_c.value"), lower: column(row, "T_c.lower"), upper: column(row, "T_c.upper") },
    ne: { value: column(row, "N_e.value"), lower: column(row, "N_e.lower"), upper: column(row, "N_e.upper") },
  }));
}

export function readClockCsv(text: string): ClockRow[] {
  return csvParse(text).map((row) => ({
    name: textColumn(row, "name"),
    div: column(row, "div"),
    date: optionalColumn(row, "date"),
    predictedDate: column(row, "predicted_date"),
    deviation: optionalColumn(row, "clock_deviation"),
    outlier: textColumn(row, "is_outlier") === "true",
  }));
}

export function readTraitsCsv(text: string): Traits {
  const rows = csvParse(text);
  const [nameColumn, attribute] = rows.columns;

  if (nameColumn === undefined || attribute === undefined) {
    throw new Error("the traits table needs a node column and a trait column");
  }

  return {
    attribute,
    states: new Map(rows.map((row) => [textColumn(row, nameColumn), textColumn(row, attribute)])),
  };
}

const zClockModel = z.object({
  clock_rate: z.number(),
  intercept: z.number(),
  stats: z.union([z.literal("fixed"), z.object({ estimated: z.object({ r_val: z.number() }) })]),
});

const zAugurNodeData = z.object({
  clock: z.object({ rate: z.number(), rate_std: z.number().optional() }).optional(),
  nodes: z.record(z.string(), z.object({ branch_length: z.number().optional() })).optional(),
});

const zGtr = z.object({ model_name: z.string(), mu: z.number() });

function column(row: DSVRowString, name: string): number {
  return parseNumber(textColumn(row, name));
}

function optionalColumn(row: DSVRowString, name: string): number | undefined {
  return parseOptionalNumber(textColumn(row, name));
}

function textColumn(row: DSVRowString, name: string): string {
  const value = row[name];

  if (value === undefined) {
    throw new Error(`the table has no column "${name}"`);
  }

  return value;
}
