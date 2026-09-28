import type { IterationEvent, LogLevel, RunEvent, TerminalEvent } from "@neherlab/app-contracts";
import { DateTime } from "luxon";

import { fromJsonFloat } from "./numbers";

export type { RunEvent };

export type LogFilter = "all" | "warnings";

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
  seq: number;
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

export type StageState = "running" | "done" | "stopped";

export interface StageSection {
  key: string;
  name: string;
  state: StageState;
  seconds: number;
  entries: readonly LogEntry[];
  defaultOpen: boolean;
}

const PRELUDE_SEQ = -1;

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
      (filter === "all" || WARNING_LEVELS.has(entry.level)) &&
      (needle === "" || entry.message.toLowerCase().includes(needle)),
  );
}

export function stageSections(progress: RunProgress): StageSection[] {
  const logs = progress.entries.filter((entry) => entry.kind === "log");
  const first = progress.stages.at(0);
  const hasPrelude = logs.some((entry) => first === undefined || entry.seq < first.seq);

  const prelude: StageSpan = {
    seq: PRELUDE_SEQ,
    name: "Start",
    startSeconds: 0,
    endSeconds: first?.startSeconds ?? (progress.terminal === undefined ? undefined : (logs.at(-1)?.seconds ?? 0)),
  };

  const spans = hasPrelude ? [prelude, ...progress.stages] : progress.stages;
  const ended = stageEndState(progress.terminal);

  return spans.map((span, index) => {
    const next = spans.at(index + 1)?.seq ?? Number.POSITIVE_INFINITY;
    const isLast = index === spans.length - 1;
    const state = span.endSeconds === undefined ? "running" : isLast ? ended : "done";

    return {
      key: String(span.seq),
      name: span.name,
      state,
      seconds: (span.endSeconds ?? span.startSeconds) - span.startSeconds,
      entries: logs.filter((entry) => entry.seq > span.seq && entry.seq < next),
      defaultOpen: state === "running",
    };
  });
}

export function openSections(sections: readonly StageSection[], overrides: ReadonlyMap<string, boolean>): string[] {
  return sections.flatMap(({ key, defaultOpen }) => ((overrides.get(key) ?? defaultOpen) ? [key] : []));
}

export function sectionOverrides(sections: readonly StageSection[], open: readonly string[]): Map<string, boolean> {
  const opened = new Set(open);

  return new Map(
    sections.flatMap((section): Array<[string, boolean]> => {
      const isOpen = opened.has(section.key);

      return isOpen === section.defaultOpen ? [] : [[section.key, isOpen]];
    }),
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
    stages: [...closeLastStage(progress.stages, seconds), { seq, name, startSeconds: seconds, endSeconds: undefined }],
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
    : [...stages.slice(0, -1), { ...last, endSeconds: seconds }];
}

function stageEndState(terminal: TerminalEvent | undefined): StageState {
  return terminal === undefined || terminal.status === "ok" ? "done" : "stopped";
}

function stageMessage(name: string, message: string): string {
  return message === "" || message === name ? name : `${name}: ${message}`;
}

function iterationPoint(event: IterationEvent): IterationPoint {
  return {
    iteration: event.iteration,
    clockRate: fromJsonFloat(event.clock_rate),
    rSquared: optionalFloat(event.r_squared),
    maxTimeChange: optionalFloat(event.max_time_change),
    rmsTimeChange: optionalFloat(event.rms_time_change),
    logLhTotal: optionalFloat(event.log_lh_total),
  };
}

function optionalFloat(value: IterationEvent["log_lh_total"]): number | undefined {
  return value === null || value === undefined ? undefined : fromJsonFloat(value);
}
