import type {
  AppCommand,
  CheckInputsRequest,
  DatasetInfo,
  InputFactsResult,
  InputKind,
  RunInput,
} from "@neherlab/app-contracts";
import * as z from "zod";

import { COMMAND_SETTINGS } from "./catalog";
import { getAt, type JsonObject, type JsonValue } from "./json";

export interface InputAssignment {
  key: string;
  value: JsonValue;
  label: string;
}

const zPath = z.string().min(1);

const DATASET_FILES: Record<InputKind, readonly string[]> = {
  tree: ["tree.nwk"],
  alignment: ["aln.fasta.xz", "aln.fasta"],
  metadata: ["metadata.tsv", "metadata.csv"],
};

export function datasetInputs(dataDir: string, dataset: DatasetInfo, command: AppCommand): InputAssignment[] {
  return COMMAND_SETTINGS[command].inputs.flatMap((input) => {
    const file = DATASET_FILES[input.kind].find((candidate) => dataset.files.includes(candidate));

    if (file === undefined) {
      return [];
    }

    const path = joinPath(dataDir, dataset.name, file);

    return [{ key: input.kind, value: input.kind === "alignment" ? [path] : path, label: `${dataset.name}/${file}` }];
  });
}

export function runInputAssignments(inputs: readonly RunInput[]): InputAssignment[] {
  const bySetting = new Map<string, RunInput[]>();

  for (const input of inputs) {
    bySetting.set(input.setting, [...(bySetting.get(input.setting) ?? []), input]);
  }

  return [...bySetting.entries()].map(([key, files]) => ({
    key,
    value: key === "alignment" ? files.map((file) => file.path) : (files[0]?.path ?? null),
    label: files.map((file) => baseName(file.path)).join(", "),
  }));
}

export function inputFactsRequest(command: AppCommand, config: JsonObject): CheckInputsRequest | null {
  const slots = new Set(COMMAND_SETTINGS[command].inputs.map((input) => input.kind));
  const tree = slots.has("tree") ? stringOrNull(getAt(config, ["tree"])) : null;
  const metadata = slots.has("metadata") ? stringOrNull(getAt(config, ["metadata"])) : null;
  const alignment = slots.has("alignment") ? stringList(getAt(config, ["alignment"])) : [];

  if (tree === null && metadata === null && alignment.length === 0) {
    return null;
  }

  return {
    tree,
    metadata,
    alignment,
    metadata_id_columns: stringList(getAt(config, ["metadata_id_columns"])),
    metadata_delimiters: stringList(getAt(config, ["metadata_delimiters"])),
    date_column: stringOrNull(getAt(config, ["date_column"])),
  };
}

export function pathList(value: JsonValue | undefined): string[] {
  if (Array.isArray(value)) {
    return value.flatMap((item) => {
      const path = stringOrNull(item);

      return path === null ? [] : [path];
    });
  }

  const path = stringOrNull(value);

  return path === null ? [] : [path];
}

export function baseName(path: string): string {
  return path.slice(Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\")) + 1);
}

function joinPath(...parts: string[]): string {
  return parts.filter((part) => part !== "").join("/");
}

function stringOrNull(value: JsonValue | undefined): string | null {
  const parsed = zPath.safeParse(value);

  return parsed.success ? parsed.data : null;
}

function stringList(value: JsonValue | undefined): string[] {
  return Array.isArray(value) ? value.map(String) : [];
}

export function slotFactsText(slot: InputKind, facts: InputFactsResult | undefined, usesDates: boolean): string | null {
  if (slot === "tree") {
    const tree = facts?.tree;

    return tree === null || tree === undefined
      ? null
      : `${tree.tips} tips, ${tree.internal_nodes} internal nodes, ${tree.polytomies} polytomies`;
  }

  if (slot === "alignment") {
    const alignment = facts?.alignment;

    if (alignment === null || alignment === undefined) {
      return null;
    }

    const sites =
      alignment.min_length === alignment.max_length
        ? `${alignment.max_length.toLocaleString("en-US")} sites`
        : `lengths ${alignment.min_length} to ${alignment.max_length}, so it is not aligned`;

    return `${alignment.sequences} sequences, ${sites}`;
  }

  const metadata = facts?.metadata;

  if (metadata === null || metadata === undefined) {
    return null;
  }

  const dateColumn =
    metadata.date_column === null || metadata.date_column === undefined
      ? "no date column found"
      : `date column ${metadata.date_column}`;

  const dates = usesDates ? `; ${dateColumn}` : "";

  return `${metadata.rows} rows; ID column ${metadata.id_column}${dates}`;
}

export function slotProblem(slot: InputKind, facts: InputFactsResult | undefined): string | null {
  return facts?.problems.find((problem) => problem.input === slot)?.message ?? null;
}
