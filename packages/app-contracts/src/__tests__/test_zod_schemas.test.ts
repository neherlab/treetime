import { describe, expect, test } from "vitest";
import { z } from "zod";

import {
  zAncestralConfig,
  zAppEvent,
  zCheckConfigResponse,
  zClockConfig,
  zCoalescentPrior,
  zCommandResults,
  zCreateRunRequest,
  zErrorResponse,
  zJobEvent,
  zLogEvent,
  zLogLevel,
  zMugrationConfig,
  zOptimizeConfig,
  zPruneConfig,
  zRunConfigResponse,
  zRunEvent,
  zSettingDifference,
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

  test.each([
    { name: "TimetreeConfig", parse: (value: unknown) => zTimetreeConfig.parse(value).seed },
    { name: "AncestralConfig", parse: (value: unknown) => zAncestralConfig.parse(value).seed },
    { name: "ClockConfig", parse: (value: unknown) => zClockConfig.parse(value).seed },
    { name: "MugrationConfig", parse: (value: unknown) => zMugrationConfig.parse(value).seed },
  ])("$name parses the 64-bit seed as a number that JSON can send", ({ parse }) => {
    const seed = parse({ seed: 42 });
    expect(JSON.stringify({ seed })).toBe('{"seed":42}');
  });

  test("timetree config rejects a seed beyond the largest exact JSON integer", () => {
    expect(zTimetreeConfig.safeParse({ seed: Number.MAX_SAFE_INTEGER + 2 }).success).toBe(false);
  });

  test("timetree config rejects a negative seed", () => {
    expect(zTimetreeConfig.safeParse({ seed: -1 }).success).toBe(false);
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

describe("zod_schemas tagged unions", () => {
  test.each([
    { name: "AppEvent", schema: zAppEvent, tag: "kind" },
    { name: "CheckConfigResponse", schema: zCheckConfigResponse, tag: "status" },
    { name: "CoalescentPrior", schema: zCoalescentPrior, tag: "kind" },
    { name: "CommandResults", schema: zCommandResults, tag: "command" },
    { name: "JobEvent", schema: zJobEvent, tag: "type" },
    { name: "RunConfigResponse", schema: zRunConfigResponse, tag: "status" },
    { name: "RunEvent", schema: zRunEvent, tag: "type" },
    { name: "SettingDifference", schema: zSettingDifference, tag: "kind" },
    { name: "TerminalEvent", schema: zTerminalEvent, tag: "status" },
  ])("$name is a discriminated union on `$tag`", ({ schema, tag }) => {
    expect(schema).toBeInstanceOf(z.ZodDiscriminatedUnion);
    expect(schema.def.discriminator).toBe(tag);
  });

  test("an unknown tag is reported at the tag with the known values", () => {
    expect(zTerminalEvent.safeParse({ status: "done", job_id: "j" }).error?.issues).toMatchObject([
      { code: "invalid_union", path: ["status"], options: ["ok", "error", "cancelled", "interrupted"] },
    ]);
  });

  test("a run event is checked against the variant its type names", () => {
    const event = { seq: 0, time: "t", type: "terminal", data: { level: "info", message: "x" } };

    expect(zRunEvent.safeParse(event).error?.issues).toMatchObject([{ path: ["data", "status"] }]);
  });

  test("a run event variant carries the sequence number and time of the event", () => {
    expect(zRunEvent.safeParse({ type: "log", data: { level: "info", message: "x" } }).success).toBe(false);
    expect(zRunEvent.safeParse({ seq: 0, time: "t", type: "log", data: { level: "info", message: "x" } }).success).toBe(
      true,
    );
  });

  test("an iteration event reads non-finite values", () => {
    const event = {
      seq: 7,
      time: "2026-09-25T10:00:00Z",
      type: "iteration",
      data: {
        iteration: 2,
        n_diff: 0,
        n_resolved: 0,
        max_time_change: 0.1,
        rms_time_change: 0.01,
        log_lh_seq: null,
        log_lh_pos: -12.5,
        log_lh_coal: "inf",
        log_lh_total: "inf",
        clock_rate: 0.001,
        r_squared: 0.8,
      },
    };

    expect(zRunEvent.parse(event)).toStrictEqual(event);
  });

  test("a variant needs the fields of its tag", () => {
    expect(zCoalescentPrior.safeParse({ kind: "fixed" }).success).toBe(false);
    expect(zCoalescentPrior.safeParse({ kind: "fixed", tc: 0.5 }).success).toBe(true);
  });

  test("an unknown field in a closed request object is rejected", () => {
    expect(zCreateRunRequest.safeParse({ command: "clock", config: {} }).success).toBe(true);
    expect(zCreateRunRequest.safeParse({ command: "clock", config: {}, extra: 1 }).error?.issues).toMatchObject([
      { code: "unrecognized_keys", keys: ["extra"] },
    ]);
  });
});

describe("zod_schemas error and enum shapes", () => {
  test("error response accepts a well-formed error", () => {
    expect(zErrorResponse.safeParse({ code: "invalid_request", message: "bad input", causes: [] }).success).toBe(true);
  });

  test("error response rejects a missing message", () => {
    expect(zErrorResponse.safeParse({ code: "invalid_request", causes: [] }).success).toBe(false);
  });

  test("error response rejects a code the back end does not send", () => {
    expect(zErrorResponse.safeParse({ code: "E_BAD", message: "bad input", causes: [] }).success).toBe(false);
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
