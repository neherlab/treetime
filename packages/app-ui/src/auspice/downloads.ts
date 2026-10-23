export function downloadPrefix(runs: FigureRuns): string {
  if (runs.length === 1) {
    const [run] = runs;

    return `treetime-${run.command}-${run.id}`;
  }

  const [first, second] = runs;

  return `treetime-compare-${first.id}-${second.id}`;
}

export type FigureRuns = readonly [FigureRun] | readonly [FigureRun, FigureRun];

export interface FigureRun {
  id: string;
  title: string;
  command: string;
}
