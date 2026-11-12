import type { RunFile } from "@neherlab/app-contracts";

export type DownloadEntry = { kind: "package" } | { kind: "file"; file: RunFile } | { kind: "all"; count: number };

export function downloadEntries(files: readonly RunFile[]): DownloadEntry[] {
  return [
    { kind: "package" },
    ...files.flatMap((file): DownloadEntry[] => (file.main ? [{ kind: "file", file }] : [])),
    { kind: "all", count: files.length },
  ];
}

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
