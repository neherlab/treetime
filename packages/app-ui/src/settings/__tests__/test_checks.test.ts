import type { CheckConfigResult, InputFactsResult } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import { configCheckMessages, formChecks, hasBlockingCheck, type CheckContext } from "../checks";
import type { InputSlotKey } from "../commands";
import { defaultConfig } from "../config";
import { setAt } from "../json";
import { commandSettings } from "../schema";

const TIMETREE_DEFAULTS = defaultConfig(commandSettings("timetree").specs);

const ALL_SLOTS: ReadonlySet<InputSlotKey> = new Set(["tree", "alignment", "metadata"]);

const ZIKA_86_FACTS: InputFactsResult = {
  tree: { tips: 86, internal_nodes: 61, polytomies: 15, unnamed_tips: 0, duplicate_tip_names: [] },
  alignment: { sequences: 86, min_length: 10807, max_length: 10807, duplicate_names: [] },
  metadata: {
    rows: 86,
    columns: ["name", "date", "country"],
    id_column: "name",
    date_column: "date",
    dates: { readable: 86, unreadable: [], exact_days: 86, on_day_1_or_15: 54 },
  },
  tips_without_metadata: [],
  tips_without_sequence: [],
  problems: [],
};

const VALID: CheckConfigResult = { status: "valid", config: {} };

function summary(overrides: Partial<CheckContext>) {
  return formChecks(context(overrides)).map((check) => [check.level, check.text]);
}

function context(overrides: Partial<CheckContext>): CheckContext {
  return {
    command: "timetree",
    config: TIMETREE_DEFAULTS,
    filledSlots: ALL_SLOTS,
    facts: ZIKA_86_FACTS,
    configCheck: VALID,
    ...overrides,
  };
}

describe("check-config presentation", () => {
  test("a valid config gives no message", () => {
    expect(configCheckMessages(VALID, false)).toStrictEqual([]);
  });

  test("each problem is one message with its help", () => {
    const result: CheckConfigResult = {
      status: "invalid",
      message: "invalid configuration: 2 problems",
      causes: [],
      problems: [
        { code: "config::enum", message: "`margnal` is not a valid value", help: "did you mean `marginal`?" },
        { code: "config::unknown-field", message: "unknown field `bogus`", help: null },
      ],
      rendered: null,
    };

    expect(configCheckMessages(result, false)).toStrictEqual([
      "`margnal` is not a valid value (did you mean `marginal`?)",
      "unknown field `bogus`",
    ]);
  });

  test("an error without problems shows the message and its causes", () => {
    const result: CheckConfigResult = {
      status: "invalid",
      message: "When resolving the arguments",
      causes: ["the attribute is required"],
      problems: [],
    };

    expect(configCheckMessages(result, false)).toStrictEqual([
      "When resolving the arguments: the attribute is required",
    ]);
  });

  test("an error without problems is left to the input checks while inputs are missing", () => {
    const result: CheckConfigResult = { status: "invalid", message: "missing tree", causes: [], problems: [] };

    expect(configCheckMessages(result, true)).toStrictEqual([]);
  });

  test("a rejected config blocks the run", () => {
    const configCheck: CheckConfigResult = {
      status: "invalid",
      message: "unknown field `x`",
      causes: [],
      problems: [],
    };

    expect(summary({ configCheck, facts: undefined })).toStrictEqual([["block", "unknown field `x`"]]);
  });
});

describe("form checks", () => {
  test("zika 86 with defaults only advises on dates rounded to the month", () => {
    expect(summary({})).toStrictEqual([
      [
        "advice",
        "54 of 86 dates fall on the 1st or 15th of a month. If only the month is known, write it as 2015-06-XX so TreeTime treats it as a range.",
      ],
    ]);
  });

  test("missing required inputs block the run", () => {
    const checks = formChecks(context({ filledSlots: new Set(["alignment"]), facts: undefined }));

    expect({ texts: checks.map((check) => check.text), blocking: hasBlockingCheck(checks) }).toStrictEqual({
      texts: ["Add a tree file.", "Add a metadata file."],
      blocking: true,
    });
  });

  test("tips without a sequence block and tips without metadata warn", () => {
    const facts: InputFactsResult = {
      ...ZIKA_86_FACTS,
      metadata: null,
      tips_without_sequence: ["a", "b", "c", "d"],
      tips_without_metadata: ["e"],
    };

    expect(summary({ facts })).toStrictEqual([
      ["block", "4 of 86 tree tips have no sequence in the alignment: a, b, c, and 1 more"],
      ["warn", "1 of 86 tree tips have no metadata row and get no date: e"],
    ]);
  });

  test("unreadable dates warn and a missing date column blocks", () => {
    const facts: InputFactsResult = {
      ...ZIKA_86_FACTS,
      metadata: {
        rows: 3,
        columns: ["name", "when"],
        id_column: "name",
        date_column: null,
        dates: { readable: 1, unreadable: ["x"], exact_days: 1, on_day_1_or_15: 0 },
      },
    };

    expect(summary({ facts })).toStrictEqual([
      ["block", "The metadata has no date column. Set the date column to one of: name, when."],
      ["warn", "1 dates cannot be read and those samples get no date: x. Use 2015-06-21, 2015-06-XX or 2015.47."],
    ]);
  });

  test("an unreadable input blocks the run", () => {
    const facts: InputFactsResult = { problems: [{ input: "tree", message: "unexpected end of file" }] };

    expect(summary({ facts })).toStrictEqual([["block", "The tree cannot be read: unexpected end of file"]]);
  });

  test("date intervals without rate uncertainty warn and offer covariation", () => {
    const config = setAt(TIMETREE_DEFAULTS, ["confidence"], true);
    const check = formChecks(context({ config, facts: undefined }))[0];

    expect({ level: check?.level, fix: check?.fix }).toStrictEqual({
      level: "warn",
      fix: { label: "Use covariation", path: ["covariation"], value: true },
    });
  });

  test("date intervals with a clock rate std. dev. give no warning", () => {
    const config = setAt(setAt(TIMETREE_DEFAULTS, ["confidence"], true), ["clock_std_dev"], 0.0001);

    expect(summary({ config, facts: undefined })).toStrictEqual([]);
  });

  test("mugration does not check dates", () => {
    const config = defaultConfig(commandSettings("mugration").specs);

    expect(summary({ command: "mugration", config })).toStrictEqual([]);
  });
});
