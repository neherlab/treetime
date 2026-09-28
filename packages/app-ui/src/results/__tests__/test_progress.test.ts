import { describe, expect, test } from "vitest";

import {
  EMPTY_PROGRESS,
  countWarnings,
  filterLog,
  foldRunEvents,
  openSections,
  sectionOverrides,
  stageSections,
  type RunEvent,
} from "../progress";

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
      { seq: 1, name: "Reading input", startSeconds: 1, endSeconds: 4 },
      { seq: 4, name: "Time tree", startSeconds: 4, endSeconds: 7 },
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

  test("warnings and search select their lines", () => {
    expect({
      warnings: filterLog(entries, "warnings", "").map((entry) => entry.seq),
      search: filterLog(entries, "all", "  SEQUENCES ").map((entry) => entry.seq),
      count: countWarnings(entries),
    }).toStrictEqual({ warnings: [5], search: [2], count: 1 });
  });
});

describe("stage sections", () => {
  const LIVE: RunEvent[] = [
    { seq: 0, time: "2026-09-25T10:00:00.000Z", type: "started", data: { job_id: "run", command: "ancestral" } },
    { seq: 1, time: "2026-09-25T10:00:00.500Z", type: "log", data: { level: "info", message: "Using 1 thread" } },
    {
      seq: 2,
      time: "2026-09-25T10:00:01.000Z",
      type: "progress",
      data: { stage: "Marginal reconstruction", fraction: 0.4, message: "" },
    },
    { seq: 3, time: "2026-09-25T10:00:02.000Z", type: "log", data: { level: "warn", message: "Short branch" } },
    {
      seq: 4,
      time: "2026-09-25T10:00:03.000Z",
      type: "progress",
      data: { stage: "Reconstructing sequences", fraction: 0.6, message: "" },
    },
    {
      seq: 5,
      time: "2026-09-25T10:00:04.000Z",
      type: "progress",
      data: { stage: "Marginal reconstruction", fraction: 0.7, message: "" },
    },
    { seq: 6, time: "2026-09-25T10:00:05.000Z", type: "log", data: { level: "info", message: "Converged" } },
  ];

  test("each stage holds the log lines emitted while it ran, and lines before the first stage form a start section", () => {
    const sections = stageSections(foldRunEvents(EMPTY_PROGRESS, LIVE));

    expect(
      sections.map((section) => ({
        key: section.key,
        name: section.name,
        state: section.state,
        seconds: section.seconds,
        lines: section.entries.map((entry) => entry.message),
      })),
    ).toStrictEqual([
      { key: "-1", name: "Start", state: "done", seconds: 1, lines: ["Using 1 thread"] },
      { key: "2", name: "Marginal reconstruction", state: "done", seconds: 2, lines: ["Short branch"] },
      { key: "4", name: "Reconstructing sequences", state: "done", seconds: 1, lines: [] },
      { key: "5", name: "Marginal reconstruction", state: "running", seconds: 0, lines: ["Converged"] },
    ]);
  });

  test("no start section when the first line follows the first stage", () => {
    const sections = stageSections(foldRunEvents(EMPTY_PROGRESS, EVENTS));

    expect(sections.map((section) => [section.name, section.entries.length])).toStrictEqual([
      ["Reading input", 1],
      ["Time tree", 1],
    ]);
  });

  test("the last stage of a run that did not finish is stopped, and earlier stages are done", () => {
    const sections = stageSections(foldRunEvents(EMPTY_PROGRESS, EVENTS));

    expect(sections.map((section) => section.state)).toStrictEqual(["done", "stopped"]);
  });

  test("every stage of a finished run is done", () => {
    const sections = stageSections(
      foldRunEvents(EMPTY_PROGRESS, [
        ...LIVE,
        {
          seq: 7,
          time: "2026-09-25T10:00:06.000Z",
          type: "terminal",
          data: { job_id: "run", status: "ok", result: { command: "ancestral", output_files: [] } },
        },
      ]),
    );

    expect(sections.map((section) => section.state)).toStrictEqual(["done", "done", "done", "done"]);
  });

  test("no stage of an ended run is open by default", () => {
    const sections = stageSections(
      foldRunEvents(EMPTY_PROGRESS, [
        ...LIVE,
        { seq: 7, time: "2026-09-25T10:00:06.000Z", type: "terminal", data: { job_id: "run", status: "cancelled" } },
      ]),
    );

    expect(openSections(sections, new Map())).toStrictEqual([]);
  });

  test("lines before any stage form a running start section", () => {
    const sections = stageSections(foldRunEvents(EMPTY_PROGRESS, LIVE.slice(0, 2)));

    expect(sections.map((section) => [section.name, section.state])).toStrictEqual([["Start", "running"]]);
  });

  test("no sections before the first line or stage", () => {
    expect(stageSections(foldRunEvents(EMPTY_PROGRESS, LIVE.slice(0, 1)))).toStrictEqual([]);
  });

  describe("open sections", () => {
    const sections = stageSections(foldRunEvents(EMPTY_PROGRESS, LIVE));

    test("only the current stage is open by default", () => {
      expect(openSections(sections, new Map())).toStrictEqual(["5"]);
    });

    test("opening an earlier stage and closing the current one are kept as overrides", () => {
      const overrides = sectionOverrides(sections, ["2"]);

      expect({ overrides: [...overrides], open: openSections(sections, overrides) }).toStrictEqual({
        overrides: [
          ["2", true],
          ["5", false],
        ],
        open: ["2"],
      });
    });

    test("a new stage opens, the previous current stage closes, and user choices on older stages stay", () => {
      const overrides = sectionOverrides(sections, ["2", "5"]);

      const next = stageSections(
        foldRunEvents(EMPTY_PROGRESS, [
          ...LIVE,
          {
            seq: 7,
            time: "2026-09-25T10:00:06.000Z",
            type: "progress",
            data: { stage: "Writing output", fraction: 0.9, message: "" },
          },
        ]),
      );

      expect(openSections(next, overrides)).toStrictEqual(["2", "7"]);
    });

    test("a stage the user opened stays open after the run ends, and the current stage closes", () => {
      const overrides = sectionOverrides(sections, ["2", "5"]);

      const ended = stageSections(
        foldRunEvents(EMPTY_PROGRESS, [
          ...LIVE,
          { seq: 7, time: "2026-09-25T10:00:06.000Z", type: "terminal", data: { job_id: "run", status: "cancelled" } },
        ]),
      );

      expect(openSections(ended, overrides)).toStrictEqual(["2"]);
    });

    test("the default open state needs no overrides", () => {
      expect([...sectionOverrides(sections, ["5"])]).toStrictEqual([]);
    });
  });
});
