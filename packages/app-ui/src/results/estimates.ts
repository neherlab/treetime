import { daysBetween } from "../format";
import { indexClades, matchAncestors } from "./clades";
import type { AugurClock, ClockModel, TraceRow } from "./readers";
import type { ResultTree } from "./tree";

const INTERVAL_EDGE_FRACTION = 0.05;

export interface TimetreeEstimates {
  rootDate: number | undefined;
  rootInterval: readonly [number, number] | undefined;
  rootNearIntervalEdge: boolean;
  rate: number | undefined;
  rateStd: number | undefined;
  rateFixed: boolean;
  r: number | undefined;
  samples: number;
  excludedSamples: number;
  logLikelihood: number | undefined;
  iterations: number;
}

export interface TimetreeOutputs {
  tree: ResultTree;
  clockModel: ClockModel | undefined;
  augurClock: AugurClock | undefined;
  trace: readonly TraceRow[] | undefined;
}

export interface EstimateComparison {
  rootShiftDays: number | undefined;
  intervalWidthDays: readonly [number | undefined, number | undefined];
  ratePercentChange: number | undefined;
}

export interface AncestorShift {
  name: string;
  tips: number;
  dateFirst: number;
  shiftDays: number;
}

export function timetreeEstimates({ tree, clockModel, augurClock, trace }: TimetreeOutputs): TimetreeEstimates {
  const rootInterval = usableInterval(tree.root.dateInterval);
  const rootDate = tree.root.date;

  return {
    rootDate,
    rootInterval,
    rootNearIntervalEdge:
      rootDate !== undefined && rootInterval !== undefined && nearIntervalEdge(rootDate, rootInterval),
    rate: clockModel?.rate ?? augurClock?.rate,
    rateStd: augurClock?.rateStd,
    rateFixed: clockModel?.fixed ?? false,
    r: clockModel?.r,
    samples: tree.tips.length,
    excludedSamples: tree.tips.filter((tip) => tip.excluded === true).length,
    logLikelihood: trace?.at(-1)?.logLhTotal,
    iterations: trace?.length ?? 0,
  };
}

function nearIntervalEdge(value: number, [lower, upper]: readonly [number, number]): boolean {
  const width = upper - lower;

  return width > 0 && Math.min(value - lower, upper - value) <= INTERVAL_EDGE_FRACTION * width;
}

export function compareEstimates(first: TimetreeEstimates, second: TimetreeEstimates): EstimateComparison {
  return {
    rootShiftDays:
      first.rootDate === undefined || second.rootDate === undefined
        ? undefined
        : daysBetween(first.rootDate, second.rootDate),
    intervalWidthDays: [intervalWidthDays(first.rootInterval), intervalWidthDays(second.rootInterval)],
    ratePercentChange:
      first.rate === undefined || second.rate === undefined || first.rate === 0
        ? undefined
        : ((second.rate - first.rate) / first.rate) * 100,
  };
}

export function ancestorShifts(first: ResultTree, second: ResultTree): AncestorShift[] {
  return matchAncestors(indexClades(first), indexClades(second)).flatMap(({ first: node, second: other }) =>
    node.date === undefined || other.date === undefined
      ? []
      : [
          {
            name: node.name,
            tips: countTips(node),
            dateFirst: node.date,
            shiftDays: daysBetween(node.date, other.date),
          },
        ],
  );
}

export function intervalWidthDays(interval: readonly [number, number] | undefined): number | undefined {
  return interval === undefined ? undefined : daysBetween(interval[0], interval[1]);
}

function usableInterval(interval: readonly [number, number] | undefined): readonly [number, number] | undefined {
  return interval !== undefined && interval[1] > interval[0] ? interval : undefined;
}

function countTips(node: ResultTree["root"]): number {
  return node.children.length === 0 ? 1 : node.children.reduce((sum, child) => sum + countTips(child), 0);
}
