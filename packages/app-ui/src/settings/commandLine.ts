import type { AppCommand } from "@neherlab/app-contracts";
import { quote } from "shlex";
import { stringify } from "yaml";

import { isChanged, settingValue } from "./config";
import { groupedSpecs } from "./groups";
import { isJsonObject, isString, setAt, type JsonObject, type JsonValue } from "./json";
import type { SettingSpec } from "./schema";

type CodeLineKind = "command" | "input" | "changed" | "output" | "comment";

export interface CodeLine {
  text: string;
  kind: CodeLineKind;
}

export const RUN_OUTPUT_DIR = "out";

const OUTPUT_ALL_KEY = "output_all";

const SCHEMA_URL_BASE = "https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas";

export function commandLineLines(command: AppCommand, specs: readonly SettingSpec[], config: JsonObject): CodeLine[] {
  const lines: CodeLine[] = [{ text: `treetime ${command}`, kind: "command" }];

  for (const spec of inputSpecs(specs, config)) {
    lines.push(flagLine(spec, settingValue(config, spec), "input"));
  }

  for (const spec of changedInOrder(specs, config)) {
    lines.push(flagLine(spec, settingValue(config, spec), "changed"));
  }

  lines.push({ text: `--output-all ${quote(outputDir(config))}`, kind: "output" });

  return lines;
}

export function commandLineText(lines: readonly CodeLine[]): string {
  const [first, ...rest] = lines.filter((line) => line.kind !== "comment");

  if (first === undefined) {
    return "";
  }

  return [first.text, ...rest.map((line) => `  ${line.text}`)].join(" \\\n");
}

export function yamlLines(command: AppCommand, specs: readonly SettingSpec[], config: JsonObject): CodeLine[] {
  const lines: CodeLine[] = [
    { text: `# yaml-language-server: $schema=${SCHEMA_URL_BASE}/input-config-${command}.schema.json`, kind: "comment" },
    { text: `# treetime ${command} --config run.yaml`, kind: "comment" },
  ];

  for (const spec of inputSpecs(specs, config)) {
    lines.push(...yamlEntry(spec.key, settingValue(config, spec), "input"));
  }

  let changed: JsonObject = {};

  for (const spec of changedInOrder(specs, config)) {
    changed = setAt(changed, spec.path, settingValue(config, spec));
  }

  for (const [key, value] of Object.entries(changed)) {
    lines.push(...yamlEntry(key, value, "changed"));
  }

  lines.push(...yamlEntry(OUTPUT_ALL_KEY, outputDir(config), "output"));

  return lines;
}

export function yamlText(lines: readonly CodeLine[]): string {
  return `${lines.map((line) => line.text).join("\n")}\n`;
}

export function flagTokens(spec: SettingSpec, value: JsonValue): string[] | null {
  const [, maxArgs] = spec.numArgs;

  if (spec.numArgs[0] === 0 && maxArgs === 0) {
    return value === true ? [spec.flag] : null;
  }

  if (value === null) {
    return null;
  }

  if (Array.isArray(value)) {
    return listTokens(spec, value);
  }

  return [spec.flag, cliValue(spec, value)];
}

function listTokens(spec: SettingSpec, values: readonly JsonValue[]): string[] | null {
  const items = values.map((value) => cliValue(spec, value));
  const [minArgs, maxArgs] = spec.numArgs;

  if (items.length === 0) {
    return null;
  }

  if (spec.valueDelimiter !== null) {
    return [spec.flag, items.join(spec.valueDelimiter)];
  }

  if (maxArgs === null || (items.length >= minArgs && items.length <= maxArgs)) {
    return [spec.flag, ...items];
  }

  if (items.length % maxArgs !== 0) {
    return null;
  }

  return Array.from({ length: items.length / maxArgs }, (_, index) => [
    spec.flag,
    ...items.slice(index * maxArgs, (index + 1) * maxArgs),
  ]).flat();
}

function cliValue(spec: SettingSpec, value: JsonValue): string {
  const text = Array.isArray(value) || isJsonObject(value) ? JSON.stringify(value) : scalarText(value);

  return spec.cliValues[text] ?? text;
}

function flagLine(spec: SettingSpec, value: JsonValue, kind: CodeLineKind): CodeLine {
  const tokens = flagTokens(spec, value);

  if (tokens === null) {
    return {
      text: `# ${spec.key} = ${JSON.stringify(value)} has no command-line form; use the YAML config`,
      kind: "comment",
    };
  }

  return { text: tokens.map((token) => quote(token)).join(" "), kind };
}

function yamlEntry(key: string, value: JsonValue, kind: CodeLineKind): CodeLine[] {
  return stringify({ [key]: value })
    .trimEnd()
    .split("\n")
    .map((text) => ({ text, kind }));
}

function inputSpecs(specs: readonly SettingSpec[], config: JsonObject): SettingSpec[] {
  return groupedSpecs(specs).flatMap(([, members]) =>
    members.filter((spec) => {
      if (spec.pathRole !== "input" && spec.pathRole !== "input-template") {
        return false;
      }

      const value = settingValue(config, spec);

      return Array.isArray(value) ? value.length > 0 : value !== null && value !== "";
    }),
  );
}

function changedInOrder(specs: readonly SettingSpec[], config: JsonObject): SettingSpec[] {
  return groupedSpecs(specs).flatMap(([, members]) => members.filter((spec) => isChanged(config, spec)));
}

function outputDir(config: JsonObject): string {
  const value = config[OUTPUT_ALL_KEY];

  return isString(value) ? value : RUN_OUTPUT_DIR;
}

function scalarText(value: string | number | boolean | null): string {
  return value === null ? "null" : `${value}`;
}
