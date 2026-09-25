import { describe, expect, test } from "vitest";

import {
  zAncestralConfig,
  zClockConfig,
  zErrorResponse,
  zLogEvent,
  zLogLevel,
  zMugrationConfig,
  zOptimizeConfig,
  zPruneConfig,
  zTerminalEvent,
  zTimetreeConfig,
  zVersionInfo,
} from "../generated/zod.gen";

interface SchemaCase {
  name: string;
  schema: { safeParse(value: unknown): { success: boolean } };
}

const CONFIG_CASES: SchemaCase[] = [
  { name: "TimetreeConfig", schema: zTimetreeConfig },
  { name: "OptimizeConfig", schema: zOptimizeConfig },
  { name: "PruneConfig", schema: zPruneConfig },
  { name: "AncestralConfig", schema: zAncestralConfig },
  { name: "ClockConfig", schema: zClockConfig },
  { name: "MugrationConfig", schema: zMugrationConfig },
];

describe("zod_schemas command configs", () => {
  test.each(CONFIG_CASES)("$name accepts an empty object because every setting has a default", ({ schema }) => {
    expect(schema.safeParse({}).success).toBe(true);
  });

  test.each(CONFIG_CASES)("$name accepts the input and output settings", ({ schema }) => {
    expect(schema.safeParse({ tree: "t.nwk", output_all: "out", output_selection: ["Nwk"] }).success).toBe(true);
  });

  test.each(CONFIG_CASES)("$name rejects an output selection outside its enum", ({ schema }) => {
    expect(schema.safeParse({ output_selection: ["NotAnOutput"] }).success).toBe(false);
  });

  test("timetree config fills the clock filter default", () => {
    const parsed = zTimetreeConfig.parse({});
    expect(parsed.clock_filter).toBe(3);
  });

  test("timetree config accepts nested settings", () => {
    expect(zTimetreeConfig.safeParse({ tree: "t.nwk", relax: [1, 0.5], max_iter: 4, seed: 7 }).success).toBe(true);
  });

  test("clock config accepts nested branch split settings", () => {
    expect(zClockConfig.safeParse({ branch_split: { method: "brent", brent_max_iters: 20 } }).success).toBe(true);
  });

  test("clock config rejects an unknown branch split method", () => {
    expect(zClockConfig.safeParse({ branch_split: { method: "newton" } }).success).toBe(false);
  });

  test("ancestral config rejects a wrong type", () => {
    expect(zAncestralConfig.safeParse({ tree: "t.nwk", dense: "yes" }).success).toBe(false);
  });

  test("ancestral config rejects an unknown method", () => {
    expect(zAncestralConfig.safeParse({ method_anc: "bogus" }).success).toBe(false);
  });

  test("mugration config accepts a fractional pseudocount", () => {
    expect(zMugrationConfig.safeParse({ attribute: "country", pc: 0.01 }).success).toBe(true);
  });

  test("optimize config rejects a negative iteration count", () => {
    expect(zOptimizeConfig.safeParse({ max_iter: -3 }).success).toBe(false);
  });

  test("prune config rejects a non-numeric branch length threshold", () => {
    expect(zPruneConfig.safeParse({ prune_short: "short" }).success).toBe(false);
  });
});

describe("zod_schemas terminal events", () => {
  test("an ok terminal event carries the outcome", () => {
    const event = {
      status: "ok",
      job_id: "j",
      result: { command: "clock", output_files: [{ path: "out/clock.nwk", kind: "nwk" }] },
    };

    expect(zTerminalEvent.safeParse(event).success).toBe(true);
  });

  test("an error terminal event needs the cause chain", () => {
    expect(zTerminalEvent.safeParse({ status: "error", job_id: "j", message: "failed" }).success).toBe(false);
  });

  test("a terminal event rejects an unknown status", () => {
    expect(zTerminalEvent.safeParse({ status: "done", job_id: "j" }).success).toBe(false);
  });
});

describe("zod_schemas error and enum shapes", () => {
  test("error response accepts a well-formed error", () => {
    expect(zErrorResponse.safeParse({ code: "E_BAD", message: "bad input" }).success).toBe(true);
  });

  test("error response rejects a missing message", () => {
    expect(zErrorResponse.safeParse({ code: "E_BAD" }).success).toBe(false);
  });

  test.each(["trace", "debug", "info", "warn", "error"])("LogLevel accepts %s, the spelling Rust sends", (level) => {
    expect(zLogLevel.safeParse(level).success).toBe(true);
  });

  test("log level rejects an unknown level", () => {
    expect(zLogLevel.safeParse("fatal").success).toBe(false);
  });

  test("log event rejects an out-of-alphabet level", () => {
    expect(zLogEvent.safeParse({ level: "verbose", message: "x" }).success).toBe(false);
  });

  test("version info rejects a non-string version", () => {
    expect(zVersionInfo.safeParse({ version: 1 }).success).toBe(false);
  });
});
