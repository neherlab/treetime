import type { AppCommand, OutputSelection } from "@neherlab/app-contracts";

import { mutedAuspiceDocument } from "./colors";
import { outputPath, type RunFileEntry } from "./files";
import {
  readAugurClock,
  readClockCsv,
  readClockModel,
  readCoalescentTsv,
  readGtr,
  readTotalBranchLength,
  readTraitsCsv,
  readTracelog,
  type AugurClock,
  type ClockModel,
  type ClockRow,
  type GtrSummary,
  type SkylineSegment,
  type TraceRow,
  type Traits,
} from "./readers";
import { parseAuspiceJson, readAuspiceTree, type AuspiceDocument, type ResultTree } from "./tree";

export interface RunResults {
  auspice: AuspiceOutput | undefined;
  clockModel: ClockModel | undefined;
  augurClock: AugurClock | undefined;
  totalBranchLength: number | undefined;
  gtr: GtrSummary | undefined;
  trace: TraceRow[] | undefined;
  skyline: SkylineSegment[] | undefined;
  clockRows: ClockRow[] | undefined;
  traits: Traits | undefined;
  problems: OutputProblem[];
}

export interface AuspiceOutput {
  document: AuspiceDocument;
  tree: ResultTree;
}

interface OutputProblem {
  path: string;
  message: string;
}

const READERS = {
  auspice: (text: string): Partial<RunResults> => {
    const tree = readAuspiceTree(parseAuspiceJson(text));

    return { auspice: { document: mutedAuspiceDocument(text, tree), tree } };
  },
  "clock-model": (text: string): Partial<RunResults> => ({ clockModel: readClockModel(text) }),
  "augur-node-data": (text: string): Partial<RunResults> => ({
    augurClock: readAugurClock(text),
    totalBranchLength: readTotalBranchLength(text),
  }),
  gtr: (text: string): Partial<RunResults> => ({ gtr: readGtr(text) }),
  tracelog: (text: string): Partial<RunResults> => ({ trace: readTracelog(text) }),
  "coalescent-tsv": (text: string): Partial<RunResults> => ({ skyline: readCoalescentTsv(text) }),
  "clock-csv": (text: string): Partial<RunResults> => ({ clockRows: readClockCsv(text) }),
  "traits-csv": (text: string): Partial<RunResults> => ({ traits: readTraitsCsv(text) }),
} satisfies Partial<Record<OutputSelection, (text: string) => Partial<RunResults>>>;

type ReadKind = keyof typeof READERS;

const RESULT_OUTPUTS: Readonly<Record<AppCommand, readonly ReadKind[]>> = {
  timetree: ["auspice", "clock-model", "augur-node-data", "tracelog", "coalescent-tsv"],
  clock: ["auspice", "clock-model", "clock-csv"],
  ancestral: ["auspice"],
  mugration: ["auspice", "traits-csv"],
  optimize: ["auspice", "augur-node-data", "gtr"],
  prune: ["auspice"],
};

const EMPTY_RESULTS: RunResults = {
  auspice: undefined,
  clockModel: undefined,
  augurClock: undefined,
  totalBranchLength: undefined,
  gtr: undefined,
  trace: undefined,
  skyline: undefined,
  clockRows: undefined,
  traits: undefined,
  problems: [],
};

export async function loadRunResults(
  command: AppCommand,
  files: readonly RunFileEntry[],
  read: (path: string) => Promise<string>,
): Promise<RunResults> {
  const parts = await Promise.all(
    RESULT_OUTPUTS[command].map(async (kind): Promise<Partial<RunResults>> => {
      const path = outputPath(files, kind);

      if (path === undefined) {
        return {};
      }

      try {
        return READERS[kind](await read(path));
      } catch (error: unknown) {
        return { problems: [{ path, message: error instanceof Error ? error.message : String(error) }] };
      }
    }),
  );

  const merged = parts.reduce<Partial<RunResults>>((all, part) => Object.assign(all, part), {});

  return { ...EMPTY_RESULTS, ...merged, problems: parts.flatMap((part) => part.problems ?? []) };
}
