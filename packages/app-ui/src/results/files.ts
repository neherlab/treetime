import type { RunFile } from "@neherlab/app-contracts";

export type DownloadEntry = { kind: "package" } | { kind: "file"; file: RunFile } | { kind: "all" };

export function downloadEntries(files: readonly RunFile[]): DownloadEntry[] {
  return [
    { kind: "package" },
    ...files.flatMap((file): DownloadEntry[] => (file.main ? [{ kind: "file", file }] : [])),
    { kind: "all" },
  ];
}

export function totalSize(files: readonly RunFile[]): number {
  return files.reduce((sum, file) => sum + file.size, 0);
}
