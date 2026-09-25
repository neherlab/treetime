import type { OutputSelection, Parsed, zRunFile } from "@neherlab/app-contracts";

export type RunFileEntry = Parsed<typeof zRunFile>;

const DESCRIPTIONS: Readonly<Partial<Record<OutputSelection, string>>> = {
  nwk: "Tree in Newick format",
  nexus: "Tree with node annotations in Nexus format",
  auspice: "Tree for Auspice and Nextstrain",
  "mat-pb": "Mutation-annotated tree (UShER protobuf)",
  "mat-json": "Mutation-annotated tree (JSON)",
  "graph-json": "Tree as a graph (JSON)",
  dot: "Tree as a graph (Graphviz)",
  "augur-node-data": "Node data for augur: dates, intervals, mutations",
  gtr: "Substitution model",
  "clock-model": "Clock rate, intercept and regression statistics",
  "confidence-tsv": "Date intervals of every node",
  "confidence-csv": "State probabilities of every node",
  "reconstructed-nuc-fasta": "Sequences of samples and ancestors",
  "reconstructed-aa-fasta": "Protein sequences of samples and ancestors",
  "traits-csv": "Inferred state of every node",
  "clock-csv": "Root-to-tip distance, date and clock prediction of every node",
  tracelog: "Convergence values of every iteration",
  "coalescent-tsv": "Coalescent time scale and effective population size",
  "coalescent-csv": "Coalescent time scale and effective population size",
  "coalescent-json": "Coalescent model and its likelihood",
};

export function outputPath(files: readonly RunFileEntry[], kind: OutputSelection): string | undefined {
  return files.find((file) => file.kind === kind)?.path;
}

export function fileDescription(file: RunFileEntry): string {
  return file.kind === null || file.kind === undefined ? "" : (DESCRIPTIONS[file.kind] ?? "");
}

export function totalSize(files: readonly RunFileEntry[]): number {
  return files.reduce((sum, file) => sum + file.size, 0);
}

export function downloadName(title: string, suffix: string): string {
  const stem = title
    .trim()
    .replaceAll(/[^\w.-]+/gu, "-")
    .replaceAll(/^-+|-+$/gu, "");

  return `${stem === "" ? "treetime-run" : stem}${suffix}`;
}
