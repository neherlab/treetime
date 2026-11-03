import { datasetSummary } from "auspice/src/components/info/datasetSummary";
import { MUTATION_COUNT_PREFIX } from "auspice/src/util/computeDistanceMetric";

import { translateAuspice } from "./i18n";
import { isDatedTree, type AuspiceState } from "./state";

export function downloadPrefix(runs: FigureRuns): string {
  if (runs.length === 1) {
    const [run] = runs;

    return `treetime-${run.command}-${run.id}`;
  }

  const [first, second] = runs;

  return `treetime-compare-${first.id}-${second.id}`;
}

export function treeFilePrefix(prefix: string, controls: Pick<AuspiceState["controls"], "distanceMeasure">): string {
  const { distanceMeasure } = controls;

  return distanceMeasure.startsWith(MUTATION_COUNT_PREFIX)
    ? `${prefix}_mutations-${distanceMeasure.slice(MUTATION_COUNT_PREFIX.length)}`
    : prefix;
}

export function figureSummary({
  metadata,
  tree,
  controls,
}: Pick<AuspiceState, "metadata" | "tree" | "controls">): string {
  return datasetSummary({
    nodes: tree.nodes ?? [],
    visibility: tree.visibility ?? [],
    mainTreeNumTips: metadata.mainTreeNumTips,
    temporalTree: isDatedTree(controls),
    t: translateAuspice,
  });
}

export type FigureRuns = readonly [FigureRun] | readonly [FigureRun, FigureRun];

export interface FigureRun {
  id: string;
  title: string;
  command: string;
}
