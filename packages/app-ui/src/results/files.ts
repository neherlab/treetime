import type { RunFile } from "@neherlab/app-contracts";

export type RunFileEntry = RunFile;

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
