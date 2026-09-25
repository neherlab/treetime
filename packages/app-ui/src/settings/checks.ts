import type { AppCommand, CheckConfigResult, InputFactsResult } from "@neherlab/app-contracts";

import { COMMAND_INFO, INPUT_SLOT_INFO, type InputSlotKey } from "./commands";
import { getAt, type JsonObject, type JsonValue } from "./json";

export type CheckLevel = "block" | "warn" | "advice";

interface CheckFix {
  label: string;
  path: string[];
  value: JsonValue;
}

export interface Check {
  id: string;
  level: CheckLevel;
  text: string;
  fix: CheckFix | null;
}

export interface CheckContext {
  command: AppCommand;
  config: JsonObject;
  filledSlots: ReadonlySet<InputSlotKey>;
  facts: InputFactsResult | undefined;
  configCheck: CheckConfigResult | undefined;
}

const NAMES_SHOWN = 3;

const MONTH_ROUNDING_SHARE = 0.25;

export function formChecks(context: CheckContext): Check[] {
  const missing = missingInputChecks(context);

  return [
    ...missing,
    ...inputProblemChecks(context),
    ...configChecks(context, missing.length > 0),
    ...sequenceChecks(context),
    ...metadataChecks(context),
    ...settingWarnings(context),
    ...dateAdvice(context),
  ];
}

export function hasBlockingCheck(checks: readonly Check[]): boolean {
  return checks.some((check) => check.level === "block");
}

export function configCheckMessages(result: CheckConfigResult, schemaProblemsOnly: boolean): string[] {
  if (result.status === "valid") {
    return [];
  }

  const problems = result.problems.map((problem) =>
    problem.help === null || problem.help === undefined || problem.help === ""
      ? problem.message
      : `${problem.message} (${problem.help})`,
  );

  if (problems.length > 0 || schemaProblemsOnly) {
    return problems;
  }

  return [[result.message, ...result.causes].join(": ")];
}

function missingInputChecks(context: CheckContext): Check[] {
  return COMMAND_INFO[context.command].slots.flatMap((slot) => {
    if (slot.need !== "required" || context.filledSlots.has(slot.key)) {
      return [];
    }

    return [block(`missing-${slot.key}`, `Add a ${INPUT_SLOT_INFO[slot.key].label.toLowerCase()} file.`)];
  });
}

function inputProblemChecks(context: CheckContext): Check[] {
  return (context.facts?.problems ?? []).map((problem) =>
    block(`unreadable-${problem.input}`, `The ${problem.input} cannot be read: ${problem.message}`),
  );
}

function configChecks(context: CheckContext, inputsMissing: boolean): Check[] {
  if (context.configCheck === undefined) {
    return [];
  }

  return configCheckMessages(context.configCheck, inputsMissing).map((message, index) =>
    block(`config-${index}`, message),
  );
}

function sequenceChecks(context: CheckContext): Check[] {
  const usesAlignment = COMMAND_INFO[context.command].slots.some((slot) => slot.key === "alignment");
  const missing = context.facts?.tips_without_sequence ?? [];
  const tips = context.facts?.tree?.tips;

  if (!usesAlignment || missing.length === 0 || tips === undefined) {
    return [];
  }

  return [
    block(
      "tips-without-sequence",
      `${missing.length} of ${tips} tree tips have no sequence in the alignment: ${nameList(missing)}`,
    ),
  ];
}

function metadataChecks(context: CheckContext): Check[] {
  const facts = context.facts;
  const usesMetadata = COMMAND_INFO[context.command].slots.some((slot) => slot.key === "metadata");

  if (facts === undefined || !usesMetadata) {
    return [];
  }

  const checks: Check[] = [];
  const missing = facts.tips_without_metadata ?? [];
  const tips = facts.tree?.tips;

  if (missing.length > 0 && tips !== undefined) {
    const consequence = COMMAND_INFO[context.command].usesDates ? "get no date" : "get no trait value";
    checks.push(
      warn(
        "tips-without-metadata",
        `${missing.length} of ${tips} tree tips have no metadata row and ${consequence}: ${nameList(missing)}`,
      ),
    );
  }

  const metadata = facts.metadata;

  if (COMMAND_INFO[context.command].usesDates && metadata !== null && metadata !== undefined) {
    if (metadata.date_column === null || metadata.date_column === undefined) {
      checks.push(
        block(
          "no-date-column",
          `The metadata has no date column. Set the date column to one of: ${metadata.columns.join(", ")}.`,
        ),
      );
    }

    const unreadable = metadata.dates?.unreadable ?? [];

    if (unreadable.length > 0) {
      checks.push(
        warn(
          "unreadable-dates",
          `${unreadable.length} dates cannot be read and those samples get no date: ${nameList(unreadable)}. Use 2015-06-21, 2015-06-XX or 2015.47.`,
        ),
      );
    }
  }

  return checks;
}

function settingWarnings(context: CheckContext): Check[] {
  if (context.command !== "timetree") {
    return [];
  }

  const confidence = getAt(context.config, ["confidence"]);
  const covariation = getAt(context.config, ["covariation"]);
  const clockStdDev = getAt(context.config, ["clock_std_dev"]);

  if (confidence !== true || covariation === true || (clockStdDev !== null && clockStdDev !== undefined)) {
    return [];
  }

  return [
    {
      id: "confidence-without-rate-uncertainty",
      level: "warn",
      text: "Date intervals need rate uncertainty: without the covariation-aware regression or a clock rate std. dev., this run writes no intervals.",
      fix: { label: "Use covariation", path: ["covariation"], value: true },
    },
  ];
}

function dateAdvice(context: CheckContext): Check[] {
  const dates = context.facts?.metadata?.dates;

  if (!COMMAND_INFO[context.command].usesDates || dates === null || dates === undefined || dates.exact_days === 0) {
    return [];
  }

  if (dates.on_day_1_or_15 <= dates.exact_days * MONTH_ROUNDING_SHARE) {
    return [];
  }

  return [
    {
      id: "dates-on-day-1-or-15",
      level: "advice",
      text: `${dates.on_day_1_or_15} of ${dates.exact_days} dates fall on the 1st or 15th of a month. If only the month is known, write it as 2015-06-XX so TreeTime treats it as a range.`,
      fix: null,
    },
  ];
}

function nameList(names: readonly string[]): string {
  const shown = names.slice(0, NAMES_SHOWN).join(", ");

  return names.length > NAMES_SHOWN ? `${shown}, and ${names.length - NAMES_SHOWN} more` : shown;
}

function block(id: string, text: string): Check {
  return { id, level: "block", text, fix: null };
}

function warn(id: string, text: string): Check {
  return { id, level: "warn", text, fix: null };
}
