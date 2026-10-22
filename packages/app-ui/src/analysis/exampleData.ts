import {
  type AppCommand,
  type DataIndex,
  type ExamplesFolder,
  type IndexedConfig,
  type IndexedDataset,
  type IndexedFile,
  type InputKind,
  type SparseConfig,
  type UiDraftSource,
  zDataIndex,
} from "@neherlab/app-contracts";

import { fuzzyFilter } from "../fuzzy";
import { commandSettings } from "../settings/catalog";
import { baseName, type InputAssignment } from "../settings/inputs";
import { type StoredFile, storedSource, storedValue } from "./fileSelection";

const DATA_FOLDER = "data";

const FETCH_TIMEOUT_MS = 300_000;

export const DATA_INDEX_FILE = "index.json";

export function exampleDataUrl(base: string, path: string): string {
  const folder = base.endsWith("/") ? base : `${base}/`;

  return `${folder}${DATA_FOLDER}/${path.split("/").map(encodeURIComponent).join("/")}`;
}

export async function fetchExampleIndex(base: string, signal?: AbortSignal): Promise<DataIndex> {
  const response = await fetchExampleData(base, DATA_INDEX_FILE, signal);

  // oxlint-disable-next-line custom/no-contract-parse -- the index is a static file of the app, which the generated client does not fetch
  return zDataIndex.parse(await response.json());
}

export async function fetchExampleFile(base: string, path: string): Promise<File> {
  const response = await fetchExampleData(base, path);

  return new File([await response.blob()], baseName(path));
}

export function webExampleIndex(index: DataIndex): ExampleIndex {
  return { origin: "web", index };
}

export function localExampleIndex({ folder, index }: ExamplesFolder): ExampleIndex {
  return { origin: "local", folder, index };
}

export function indexedDatasetInputs(dataset: IndexedDataset, command: AppCommand): ExampleInput[] {
  return commandSettings(command).inputs.flatMap((slot) => {
    const input = dataset.inputs.find((candidate) => candidate.kind === slot.kind);

    return input === undefined ? [] : [{ key: slot.kind, list: slot.list, path: input.path }];
  });
}

export function localInputAssignment(folder: string, input: ExampleInput): InputAssignment {
  const path = localExamplePath(folder, input.path);

  return { key: input.key, value: input.list ? [path] : path, label: input.path };
}

export function localExamplePath(folder: string, path: string): string {
  const separator = folder.includes("\\") && !folder.includes("/") ? "\\" : "/";
  const trimmed = folder.endsWith(separator) ? folder.slice(0, -separator.length) : folder;

  return [trimmed, ...path.split("/").filter((part) => part !== "")].join(separator);
}

export function configInputs(config: IndexedConfig): ConfigInput[] {
  const specs = commandSettings(config.command).settings;

  return Object.entries(config.inputs).map(([key, paths]) => ({
    key,
    list: specs.some((spec) => spec.key === key && spec.kind === "list"),
    paths,
  }));
}

export function uploadedConfigInputs(
  inputs: ReadonlyArray<{ input: ConfigInput; stored: readonly StoredFile[] }>,
  runId: string,
): UploadedInputs {
  const uploaded = inputs.flatMap(({ input, stored }) => {
    const value = storedValue(stored, input.list);
    const source: UiDraftSource = { ...storedSource(stored), origin: "upload", run_id: runId };

    return value === undefined ? [] : [{ key: input.key, value, source }];
  });

  return {
    replace_inputs: Object.fromEntries(uploaded.map(({ key, value }) => [key, value])),
    sources: Object.fromEntries(uploaded.map(({ key, source }) => [key, source])),
  };
}

export function matchingExampleFiles(files: readonly IndexedFile[], query: string): IndexedFile[] {
  return fuzzyFilter(files, query, (file) => file.path);
}

export type ExampleIndex =
  | { origin: "web"; index: DataIndex }
  | { origin: "local"; folder: string; index: DataIndex | undefined };

export interface ExampleInput {
  key: InputKind;
  list: boolean;
  path: string;
}

export interface ConfigInput {
  key: string;
  list: boolean;
  paths: string[];
}

export interface UploadedInputs {
  replace_inputs: SparseConfig;
  sources: Record<string, UiDraftSource>;
}

async function fetchExampleData(base: string, path: string, signal?: AbortSignal): Promise<Response> {
  const timeout = AbortSignal.timeout(FETCH_TIMEOUT_MS);

  const response = await fetch(exampleDataUrl(base, path), {
    signal: signal === undefined ? timeout : AbortSignal.any([signal, timeout]),
  });

  if (!response.ok) {
    throw new Error(`The example file ${path} cannot be read: HTTP ${response.status} ${response.statusText}`);
  }

  return response;
}
