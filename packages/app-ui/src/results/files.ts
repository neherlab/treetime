import type { RunFile } from "@neherlab/app-contracts";

export function totalSize(files: readonly RunFile[]): number {
  return files.reduce((sum, file) => sum + file.size, 0);
}

export function downloadName(title: string, suffix: string): string {
  const stem = title
    .trim()
    .replaceAll(/[^\w.-]+/gu, "-")
    .replaceAll(/^-+|-+$/gu, "");

  return `${stem === "" ? "treetime-run" : stem}${suffix}`;
}
