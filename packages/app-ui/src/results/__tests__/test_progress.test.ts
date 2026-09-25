import { describe, expect, test } from "vitest";

import { EMPTY_PROGRESS, countWarnings, filterLog, foldRunEvents, type RunEvent } from "../progress";

const EVENTS: RunEvent[] = [
  { seq: 0, time: "2026-09-25T10:00:00.000Z", type: "started", data: { job_id: "run", command: "timetree" } },
  {
    seq: 1,
    time: "2026-09-25T10:00:01.000Z",
    type: "progress",
    data: { stage: "Reading input", fraction: 0.1, message: "" },
  },
  { seq: 2, time: "2026-09-25T10:00:01.500Z", type: "log", data: { level: "info", message: "Read 86 sequences" } },
  {
    seq: 3,
    time: "2026-09-25T10:00:02.000Z",
    type: "progress",
    data: { stage: "Reading input", fraction: 0.2, message: "tree" },
  },
  {
    seq: 4,
    time: "2026-09-25T10:00:04.000Z",
    type: "progress",
    data: { stage: "Time tree", fraction: 0.5, message: "" },
  },
  {
    seq: 5,
    time: "2026-09-25T10:00:05.000Z",
    type: "log",
    data: { level: "warn", message: "Coalescent likelihood is not finite" },
  },
  {
    seq: 6,
    time: "2026-09-25T10:00:06.000Z",
    type: "iteration",
    data: {
      iteration: 0,
      n_diff: 1,
      n_resolved: 0,
      max_time_change: 0.5,
      rms_time_change: 0.1,
      log_lh_seq: -100,
      log_lh_pos: "-inf",
      log_lh_coal: "nan",
      log_lh_total: "inf",
      clock_rate: 0.001,
      r_squared: null,
    },
  },
  { seq: 7, time: "2026-09-25T10:00:07.000Z", type: "terminal", data: { job_id: "run", status: "cancelled" } },
];

describe("run event folding", () => {
  test("stages keep the time they started and ended", () => {
    const progress = foldRunEvents(EMPTY_PROGRESS, EVENTS);

    expect(progress.stages).toStrictEqual([
      { name: "Reading input", startSeconds: 1, endSeconds: 4 },
      { name: "Time tree", startSeconds: 4, endSeconds: 7 },
    ]);
  });

  test("non-finite iteration values become infinities and NaN, missing values stay absent", () => {
    const [point] = foldRunEvents(EMPTY_PROGRESS, EVENTS).iterations;

    expect(point).toStrictEqual({
      iteration: 0,
      clockRate: 0.001,
      rSquared: undefined,
      maxTimeChange: 0.5,
      rmsTimeChange: 0.1,
      logLhTotal: Number.POSITIVE_INFINITY,
    });
  });

  test("replaying events from an earlier offset does not duplicate them", () => {
    const once = foldRunEvents(EMPTY_PROGRESS, EVENTS);
    const twice = foldRunEvents(once, EVENTS);

    expect(twice).toStrictEqual(once);
  });

  test("the fraction never decreases and the terminal event is kept", () => {
    const progress = foldRunEvents(EMPTY_PROGRESS, [
      ...EVENTS.slice(0, 5),
      {
        seq: 5,
        time: "2026-09-25T10:00:05.000Z",
        type: "progress",
        data: { stage: "Time tree", fraction: 0.3, message: "" },
      },
      ...EVENTS.slice(7),
    ]);

    expect({ fraction: progress.fraction, terminal: progress.terminal?.status }).toStrictEqual({
      fraction: 0.5,
      terminal: "cancelled",
    });
  });
});

describe("log filters", () => {
  const entries = foldRunEvents(EMPTY_PROGRESS, EVENTS).entries;

  test("the log holds stage changes and log lines in event order", () => {
    expect(entries.map((entry) => [entry.kind, entry.message])).toStrictEqual([
      ["stage", "Reading input"],
      ["log", "Read 86 sequences"],
      ["stage", "Time tree"],
      ["log", "Coalescent likelihood is not finite"],
    ]);
  });

  test("warnings, stages and search select their lines", () => {
    expect({
      warnings: filterLog(entries, "warnings", "").map((entry) => entry.seq),
      stages: filterLog(entries, "stages", "").map((entry) => entry.seq),
      search: filterLog(entries, "all", "  SEQUENCES ").map((entry) => entry.seq),
      count: countWarnings(entries),
    }).toStrictEqual({ warnings: [5], stages: [1, 4], search: [2], count: 1 });
  });
});
