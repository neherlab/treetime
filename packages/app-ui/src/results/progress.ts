import type { Parsed, zIterationEvent, zLogLevel, zRunEvent, zTerminalEvent } from "@neherlab/app-contracts";
import { DateTime } from "luxon";

import { fromJsonFloat } from "./numbers";

export type RunEvent = Parsed<typeof zRunEvent>;

type TerminalEvent = Parsed<typeof zTerminalEvent>;

type LogLevel = Parsed<typeof zLogLevel>;

export type LogFilter = "all" | "warnings" | "stages";

export interface RunProgress {
  startedAt: DateTime | undefined;
  stages: readonly StageSpan[];
  fraction: number;
  message: string;
  iterations: readonly IterationPoint[];
  entries: readonly LogEntry[];
  terminal: TerminalEvent | undefined;
  nextSeq: number;
}

interface StageSpan {
  name: string;
  startSeconds: number;
  endSeconds: number | undefined;
}

export interface IterationPoint {
  iteration: number;
  clockRate: number;
  rSquared: number | undefined;
  maxTimeChange: number | undefined;
  rmsTimeChange: number | undefined;
  logLhTotal: number | undefined;
}

export interface LogEntry {
  seq: number;
  seconds: number;
  kind: "log" | "stage";
  level: LogLevel;
  message: string;
}

export const EMPTY_PROGRESS: RunProgress = {
  startedAt: undefined,
  stages: [],
  fraction: 0,
  message: "",
  iterations: [],
  entries: [],
  terminal: undefined,
  nextSeq: 0,
};

const WARNING_LEVELS = new Set<LogLevel>(["warn", "error"]);

export function foldRunEvents(progress: RunProgress, events: readonly RunEvent[]): RunProgress {
  return events.reduce(foldRunEvent, progress);
}

export function filterLog(entries: readonly LogEntry[], filter: LogFilter, query: string): LogEntry[] {
  const needle = query.trim().toLowerCase();

  return entries.filter(
    (entry) =>
      (filter === "all" ||
        (filter === "warnings" && WARNING_LEVELS.has(entry.level)) ||
        (filter === "stages" && entry.kind === "stage")) &&
      (needle === "" || entry.message.toLowerCase().includes(needle)),
  );
}

export function countWarnings(entries: readonly LogEntry[]): number {
  return entries.filter((entry) => WARNING_LEVELS.has(entry.level)).length;
}

export function logText(entries: readonly LogEntry[]): string {
  return entries
    .map((entry) => `${formatSeconds(entry.seconds)} ${entry.level.toUpperCase()} ${entry.message}`)
    .join("\n");
}

export function formatSeconds(seconds: number): string {
  return `${seconds.toFixed(1).padStart(6)} s`;
}

function foldRunEvent(progress: RunProgress, event: RunEvent): RunProgress {
  if (event.seq < progress.nextSeq) {
    return progress;
  }

  const time = DateTime.fromISO(event.time);
  const startedAt = progress.startedAt ?? time;
  const seconds = Math.max(0, time.diff(startedAt, "seconds").seconds);
  const next = { ...progress, startedAt, nextSeq: event.seq + 1 };

  if (event.type === "progress") {
    return foldStage(next, event.seq, seconds, event.data.stage, event.data.fraction, event.data.message);
  }

  if (event.type === "log") {
    const entry: LogEntry = {
      seq: event.seq,
      seconds,
      kind: "log",
      level: event.data.level,
      message: event.data.message,
    };

    return { ...next, entries: [...next.entries, entry] };
  }

  if (event.type === "iteration") {
    return { ...next, iterations: [...next.iterations, iterationPoint(event.data)] };
  }

  if (event.type === "terminal") {
    return { ...next, terminal: event.data, stages: closeLastStage(next.stages, seconds) };
  }

  return next;
}

function foldStage(
  progress: RunProgress,
  seq: number,
  seconds: number,
  name: string,
  fraction: number,
  message: string,
): RunProgress {
  const current = progress.stages.at(-1);
  const updated = { ...progress, fraction: Math.max(progress.fraction, fraction), message };

  if (current?.name === name) {
    return updated;
  }

  return {
    ...updated,
    stages: [...closeLastStage(progress.stages, seconds), { name, startSeconds: seconds, endSeconds: undefined }],
    entries: [
      ...progress.entries,
      { seq, seconds, kind: "stage", level: "info", message: stageMessage(name, message) },
    ],
  };
}

function closeLastStage(stages: readonly StageSpan[], seconds: number): StageSpan[] {
  const last = stages.at(-1);

  return last === undefined || last.endSeconds !== undefined
    ? [...stages]
    : [...stages.slice(0, -1), { name: last.name, startSeconds: last.startSeconds, endSeconds: seconds }];
}

function stageMessage(name: string, message: string): string {
  return message === "" || message === name ? name : `${name}: ${message}`;
}

function iterationPoint(event: Parsed<typeof zIterationEvent>): IterationPoint {
  return {
    iteration: event.iteration,
    clockRate: fromJsonFloat(event.clock_rate),
    rSquared: optionalFloat(event.r_squared),
    maxTimeChange: optionalFloat(event.max_time_change),
    rmsTimeChange: optionalFloat(event.rms_time_change),
    logLhTotal: optionalFloat(event.log_lh_total),
  };
}

function optionalFloat(value: Parsed<typeof zIterationEvent>["log_lh_total"]): number | undefined {
  return value === null || value === undefined ? undefined : fromJsonFloat(value);
}
