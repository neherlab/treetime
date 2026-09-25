import type {
  AppCommand,
  CheckInputsRequest,
  Dataset,
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

export function datasetInputs(dataset: Dataset, command: AppCommand): InputAssignment[] {
  return COMMAND_SETTINGS[command].inputs.flatMap((slot) => {
    const input = dataset.inputs.find((candidate) => candidate.kind === slot.kind);

    return input === undefined
      ? []
      : [{ key: slot.kind, value: slot.list ? [input.path] : input.path, label: input.file }];
  });
}

export function runInputAssignments(command: AppCommand, inputs: readonly RunInput[]): InputAssignment[] {
  const bySetting = new Map<string, RunInput[]>();

  for (const input of inputs) {
    bySetting.set(input.setting, [...(bySetting.get(input.setting) ?? []), input]);
  }

  return [...bySetting.entries()].map(([key, files]) => ({
    key,
    value: isListInput(command, key) ? files.map((file) => file.path) : (files[0]?.path ?? null),
    label: files.map((file) => baseName(file.path)).join(", "),
  }));
}

export function inputFactsRequest(command: AppCommand, config: JsonObject): CheckInputsRequest | null {
  const given = COMMAND_SETTINGS[command].inputs.some((slot) => pathList(getAt(config, [slot.kind])).length > 0);

  return given ? { command, config } : null;
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

function isListInput(command: AppCommand, key: string): boolean {
  return COMMAND_SETTINGS[command].inputs.some((slot) => slot.kind === key && slot.list);
}

function stringOrNull(value: JsonValue | undefined): string | null {
  const parsed = zPath.safeParse(value);

  return parsed.success ? parsed.data : null;
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
