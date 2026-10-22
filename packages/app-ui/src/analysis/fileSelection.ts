import type { JsonValue, UiDraft, UiDraftSource } from "@neherlab/app-contracts";

import { pathList } from "../settings/inputs";

export interface StoredFile {
  name: string;
  path: string;
  size: number;
}

export function filesForSetting<T>(files: readonly T[], list: boolean): T[] {
  return list ? [...files] : files.slice(0, 1);
}

export async function storeInOrder<T>(
  files: readonly T[],
  store: (file: T) => Promise<StoredFile>,
  discard: (stored: readonly StoredFile[]) => void,
): Promise<StoredFile[]> {
  const stored: StoredFile[] = [];

  try {
    for (const file of files) {
      stored.push(await store(file));
    }
  } catch (error) {
    discard(stored);
    throw error;
  }

  return stored;
}

export function releasedPaths(previous: JsonValue | undefined, next: JsonValue | undefined): string[] {
  const kept = new Set(pathList(next));

  return pathList(previous).filter((path) => !kept.has(path));
}

export function draftUploadPaths({ config, sources, upload_run_id }: UiDraft, keys?: ReadonlySet<string>): string[] {
  return Object.entries(sources).flatMap(([key, source]) =>
    upload_run_id !== undefined &&
    source.origin === "upload" &&
    source.run_id === upload_run_id &&
    (keys === undefined || keys.has(key))
      ? pathList(config[key])
      : [],
  );
}

export function releasedUploads(previous: UiDraft, next: Pick<UiDraft, "config" | "sources">): string[] {
  const kept = new Set(Object.keys(next.sources).flatMap((key) => pathList(next.config[key])));

  return draftUploadPaths(previous).filter((path) => !kept.has(path));
}

export function storedValue(stored: readonly StoredFile[], list: boolean): JsonValue | undefined {
  const paths = stored.map((file) => file.path);

  if (paths.length === 0) {
    return undefined;
  }

  return list ? paths : paths[0];
}

export function storedSource(stored: readonly StoredFile[]): Pick<UiDraftSource, "label" | "size"> {
  return {
    label: stored.map((file) => file.name).join(", "),
    size: stored.reduce((total, file) => total + file.size, 0),
  };
}
