import { describe, expect, test } from "vitest";

import { APP_COMMANDS } from "../commands";
import { commandSettings, type SettingSpec } from "../schema";

function spec(command: Parameters<typeof commandSettings>[0], key: string): SettingSpec {
  const found = commandSettings(command).specs.find((candidate) => candidate.key === key);

  if (found === undefined) {
    throw new Error(`no setting ${key} in ${command}`);
  }

  return found;
}

function summary(found: SettingSpec) {
  return {
    kind: found.kind,
    nullable: found.nullable,
    defaultValue: found.defaultValue,
    flag: found.flag,
    numArgs: found.numArgs,
    valueDelimiter: found.valueDelimiter,
    pathRole: found.pathRole,
  };
}

describe("schema settings", () => {
  test.each(APP_COMMANDS)("every setting of %s has a field kind and a flag", (command) => {
    const specs = commandSettings(command).specs;
    const unannotated = specs.filter((found) => !found.flag.startsWith("--"));

    expect({ empty: specs.length === 0, unannotated }).toStrictEqual({ empty: false, unannotated: [] });
  });

  test.each([
    ["timetree", "confidence", "switch", false, false, "--confidence", [0, 0], null, null],
    ["ancestral", "dense", "tristate", true, null, "--dense", [1, 1], null, null],
    ["timetree", "max_iter", "integer", false, 2, "--max-iter", [1, 1], null, null],
    ["timetree", "seed", "integer", true, null, "--seed", [1, 1], null, null],
    ["timetree", "clock_filter", "number", false, 3, "--clock-filter", [1, 1], null, null],
    ["timetree", "clock_rate", "number", true, null, "--clock-rate", [1, 1], null, null],
    ["timetree", "date_format", "text", false, "%Y-%m-%d", "--date-format", [1, 1], null, null],
    ["timetree", "date_column", "text", true, null, "--date-column", [1, 1], null, null],
    ["timetree", "relax", "list", false, [], "--relax", [2, 2], null, null],
    ["timetree", "model_params", "list", false, [], "--model-params", [1, 1], null, null],
    ["timetree", "output_selection", "enum-list", false, [], "--output-selection", [1, 1], ",", null],
    ["timetree", "reroot", "enum", true, null, "--reroot", [1, 1], null, null],
    ["timetree", "model", "enum", false, "infer", "--model", [1, 1], null, null],
    ["optimize", "reroot", "enum", true, null, "--reroot", [0, 1], null, null],
    ["timetree", "tree", "text", true, null, "--tree", [1, 1], null, "input"],
    ["timetree", "alignment", "list", false, [], "--alignment", [1, 1], null, "input"],
    ["ancestral", "translations", "text", true, null, "--translations", [1, 1], null, "input-template"],
    ["clock", "output_all", "text", true, null, "--output-all", [1, 1], null, "output"],
    ["clock", "branch_split.n_points", "integer", false, 11, "--branch-split-grid-n-points", [1, 1], null, null],
    ["clock", "branch_split.method", "enum", false, "grid", "--branch-split-method", [1, 1], null, null],
    ["clock", "clock_regression.variance_factor", "number", false, 0, "--variance-factor", [1, 1], null, null],
  ] as const)(
    "%s %s is a %s field",
    (command, key, kind, nullable, defaultValue, flag, numArgs, valueDelimiter, pathRole) => {
      expect(summary(spec(command, key))).toStrictEqual({
        kind,
        nullable,
        defaultValue,
        flag,
        numArgs,
        valueDelimiter,
        pathRole,
      });
    },
  );

  test("nested settings keep their parent key as path", () => {
    expect(spec("clock", "branch_split.n_points").path).toStrictEqual(["branch_split", "n_points"]);
  });

  test("the model offers every named model and the inferred one with its help", () => {
    const model = spec("timetree", "model");

    expect(model.options.map((option) => option.value)).toStrictEqual([
      "jc69",
      "k80",
      "f81",
      "hky85",
      "t92",
      "tn93",
      "jtt92",
      "infer",
    ]);
    expect(model.options.at(-1)?.help).toStrictEqual(
      "Infer GTR parameters from data via Fitch parsimony substitution counts.",
    );
  });

  test("output selections map to their command-line spelling", () => {
    expect(spec("timetree", "output_selection").cliValues["AugurNodeData"]).toStrictEqual("augur-node-data");
  });

  test("a list of numbers is read as numbers", () => {
    expect(spec("timetree", "relax").itemKind).toStrictEqual("number");
  });

  test("the help is the first paragraph with its line breaks joined, and the rest is more", () => {
    const confidence = spec("timetree", "confidence");

    expect({ help: confidence.help, more: confidence.more.startsWith("`--time-marginal=always`") }).toStrictEqual({
      help: "Add rate-uncertainty to confidence intervals.",
      more: true,
    });
  });

  test("a paragraph wrapped over several lines reads as one line", () => {
    expect(spec("timetree", "confidence").more.includes("\n")).toStrictEqual(false);
  });
});
