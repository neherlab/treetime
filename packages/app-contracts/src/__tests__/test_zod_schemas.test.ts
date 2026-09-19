import { describe, expect, test } from "vitest";

import {
  zAncestralArgs,
  zClockArgs,
  zErrorResponse,
  zLogEvent,
  zLogLevel,
  zMugrationArgs,
  zOptimizeArgs,
  zPruneArgs,
  zTimetreeArgs,
  zVersionInfo,
} from "../generated/zod.gen";

interface SchemaCase {
  name: string;
  schema: { safeParse(value: unknown): { success: boolean } };
  base: Record<string, unknown>;
}

const SCHEMA_CASES: SchemaCase[] = [
  { name: "AncestralArgs", schema: zAncestralArgs, base: { tree: "t.nwk", outdir: "out" } },
  { name: "ClockArgs", schema: zClockArgs, base: { dates: "d.tsv", outdir: "out" } },
  { name: "TimetreeArgs", schema: zTimetreeArgs, base: { outdir: "out" } },
  {
    name: "MugrationArgs",
    schema: zMugrationArgs,
    base: { attribute: "country", states: "s.tsv", outdir: "out" },
  },
  { name: "OptimizeArgs", schema: zOptimizeArgs, base: { tree: "t.nwk", outdir: "out" } },
  { name: "PruneArgs", schema: zPruneArgs, base: { tree: "t.nwk", outdir: "out" } },
];

describe("zod_schemas optional fields", () => {
  test.each(SCHEMA_CASES)("$name accepts only its required fields", ({ schema, base }) => {
    expect(schema.safeParse(base).success).toBe(true);
  });

  test.each(SCHEMA_CASES)("$name accepts an empty object because every field carries a default", ({ schema }) => {
    expect(schema.safeParse({}).success).toBe(true);
  });

  test("AncestralArgs accepts optional fields when present", () => {
    const parsed = zAncestralArgs.safeParse({
      tree: "t.nwk",
      outdir: "out",
      input_fastas: ["a.fasta"],
      dense: true,
      seed: 7,
    });
    expect(parsed.success).toBe(true);
  });

  test("AncestralArgs rejects a wrong type in an optional field", () => {
    expect(zAncestralArgs.safeParse({ tree: "t.nwk", outdir: "out", dense: "yes" }).success).toBe(false);
  });
});

interface IntCase {
  name: string;
  schema: { safeParse(value: unknown): { success: boolean } };
  base: Record<string, unknown>;
  field: string;
  int64?: boolean;
}

const INT_CASES: IntCase[] = [
  { name: "AncestralArgs.gtr_iterations", schema: zAncestralArgs, base: { tree: "t", outdir: "o" }, field: "gtr_iterations" },
  { name: "AncestralArgs.seed", schema: zAncestralArgs, base: { tree: "t", outdir: "o" }, field: "seed", int64: true },
  { name: "ClockArgs.sequence_length", schema: zClockArgs, base: { dates: "d", outdir: "o" }, field: "sequence_length" },
  { name: "TimetreeArgs.max_iter", schema: zTimetreeArgs, base: { outdir: "o" }, field: "max_iter" },
  {
    name: "MugrationArgs.iterations",
    schema: zMugrationArgs,
    base: { attribute: "a", states: "s", outdir: "o" },
    field: "iterations",
  },
  { name: "OptimizeArgs.max_iter", schema: zOptimizeArgs, base: { tree: "t", outdir: "o" }, field: "max_iter" },
];

describe("zod_schemas integer fields", () => {
  test.each(INT_CASES)("$name accepts a whole number", ({ schema, base, field }) => {
    expect(schema.safeParse({ ...base, [field]: 4 }).success).toBe(true);
  });

  test.each(INT_CASES)("$name rejects a negative whole number because the Rust type is unsigned", ({ schema, base, field }) => {
    expect(schema.safeParse({ ...base, [field]: -3 }).success).toBe(false);
  });

  test.each(INT_CASES)("$name rejects a fractional number", ({ schema, base, field }) => {
    expect(schema.safeParse({ ...base, [field]: 2.5 }).success).toBe(false);
  });

  test.each(INT_CASES)("$name rejects NaN", ({ schema, base, field }) => {
    expect(schema.safeParse({ ...base, [field]: Number.NaN }).success).toBe(false);
  });

  test.each(INT_CASES.filter((c) => !c.int64))("$name rejects a numeric string", ({ schema, base, field }) => {
    expect(schema.safeParse({ ...base, [field]: "4" }).success).toBe(false);
  });

  test.each(INT_CASES.filter((c) => c.int64))(
    "$name accepts a numeric string because JSON cannot carry a 64-bit integer",
    ({ schema, base, field }) => {
      expect(schema.safeParse({ ...base, [field]: "4" }).success).toBe(true);
    },
  );
});

describe("zod_schemas number fields keep fractional values", () => {
  test("ClockArgs.clock_filter accepts a fractional number", () => {
    expect(zClockArgs.safeParse({ dates: "d", outdir: "o", clock_filter: 2.5 }).success).toBe(true);
  });

  test("MugrationArgs.pc accepts a fractional number", () => {
    expect(zMugrationArgs.safeParse({ attribute: "a", states: "s", outdir: "o", pc: 0.01 }).success).toBe(true);
  });
});

describe("zod_schemas error and enum shapes", () => {
  test("ErrorResponse accepts a well-formed error", () => {
    expect(zErrorResponse.safeParse({ code: "E_BAD", message: "bad input" }).success).toBe(true);
  });

  test("ErrorResponse rejects a missing message", () => {
    expect(zErrorResponse.safeParse({ code: "E_BAD" }).success).toBe(false);
  });

  test.each(["Trace", "Debug", "Info", "Warn", "Error"])("LogLevel accepts %s", (level) => {
    expect(zLogLevel.safeParse(level).success).toBe(true);
  });

  test("LogLevel rejects an unknown level", () => {
    expect(zLogLevel.safeParse("Fatal").success).toBe(false);
  });

  test("LogEvent rejects an out-of-alphabet level", () => {
    expect(zLogEvent.safeParse({ level: "Verbose", message: "x" }).success).toBe(false);
  });

  test("VersionInfo rejects a non-string version", () => {
    expect(zVersionInfo.safeParse({ version: 1 }).success).toBe(false);
  });
});
