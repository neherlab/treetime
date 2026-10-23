import type {
  AlignmentFacts,
  AppCommand,
  CheckInputsRequest,
  InputFacts,
  InputKind,
  JsonValue,
  RunInput,
  SparseConfig,
} from "@neherlab/app-contracts";
import * as z from "zod";

import { commandSettings } from "./catalog";

export interface InputAssignment {
  key: string;
  value: JsonValue;
  label: string;
}

const zPath = z.string().min(1);

export function runInputAssignments(command: AppCommand, inputs: readonly RunInput[]): InputAssignment[] {
  const bySetting = new Map<string, RunInput[]>();

  for (const input of inputs) {
    bySetting.set(input.setting, [...(bySetting.get(input.setting) ?? []), input]);
  }

  return [...bySetting.entries()].flatMap(([key, files]) => {
    const [first] = files;

    return first === undefined
      ? []
      : [
          {
            key,
            value: isListInput(command, key) ? files.map((file) => file.path) : first.path,
            label: files.map((file) => baseName(file.path)).join(", "),
          },
        ];
  });
}

export function inputSettingKeys(command: AppCommand): string[] {
  return commandSettings(command).settings.flatMap((spec) => (spec.role === "input" ? [spec.key] : []));
}

export function inputFactsRequest(command: AppCommand, config: SparseConfig): CheckInputsRequest | null {
  const given = commandSettings(command).inputs.some((slot) => pathList(config[slot.kind]).length > 0);

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
  return path.slice(lastSeparator(path) + 1);
}

export function folderName(path: string): string | undefined {
  const separator = lastSeparator(path);

  return separator < 0 ? undefined : path.slice(0, Math.max(separator, 1));
}

function lastSeparator(path: string): number {
  return Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\"));
}

function isListInput(command: AppCommand, key: string): boolean {
  return commandSettings(command).inputs.some((slot) => slot.kind === key && slot.list);
}

function stringOrNull(value: JsonValue | undefined): string | null {
  const parsed = zPath.safeParse(value);

  return parsed.success ? parsed.data : null;
}

export function slotFactsText(slot: InputKind, facts: InputFacts | undefined, usesDates: boolean): string[] | null {
  if (slot === "tree") {
    const tree = facts?.tree;

    return tree === undefined
      ? null
      : [`${tree.tips} tips, ${tree.internal_nodes} internal nodes, ${tree.polytomies} polytomies`];
  }

  if (slot === "alignment") {
    const alignments = facts?.alignments ?? [];

    if (alignments.length === 0) {
      return null;
    }

    const several =
      alignments.length > 1 || facts?.problems.some((problem) => problem.alignment !== undefined) === true;

    return alignments.map((alignment) => {
      const text = alignmentFactsText(alignment);

      return several ? `${alignmentName(alignment.number, alignment.path)}: ${text}` : text;
    });
  }

  const metadata = facts?.metadata;

  if (metadata === undefined) {
    return null;
  }

  const dateColumn =
    metadata.date_column === undefined ? "no date column found" : `date column ${metadata.date_column}`;

  const dates = usesDates ? `; ${dateColumn}` : "";

  return [`${metadata.rows} rows; ID column ${metadata.id_column}${dates}`];
}

export function slotProblems(slot: InputKind, facts: InputFacts | undefined): string[] {
  return (facts?.problems ?? [])
    .filter((problem) => problem.input === slot)
    .map((problem) =>
      problem.alignment === undefined ? problem.message : `Alignment ${problem.alignment}: ${problem.message}`,
    );
}

function alignmentFactsText(alignment: AlignmentFacts): string {
  const sites =
    alignment.min_length === alignment.max_length
      ? `${alignment.max_length.toLocaleString("en-US")} sites`
      : `lengths ${alignment.min_length} to ${alignment.max_length}, so it is not aligned`;

  return `${alignment.sequences} sequences, ${sites}`;
}

function alignmentName(number: number, path: string): string {
  return `Alignment ${number} (${baseName(path)})`;
}
