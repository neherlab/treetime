import { describe, expect, test } from "vitest";

import { defaultConfig } from "../config";
import { configTextCommand, withMissingInputs } from "../configText";
import { setAt } from "../json";
import { commandSettings } from "../schema";

describe("config text command", () => {
  test("the schema directive names the command", () => {
    const text =
      "# yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-ancestral.schema.json\ntree: t.nwk\n";

    expect(configTextCommand(text)).toStrictEqual("ancestral");
  });

  test("a pipeline schema names no app command", () => {
    expect(configTextCommand("# $schema=.../input-config-pipeline.schema.json\n")).toStrictEqual(null);
  });

  test("a config without a directive names no command", () => {
    expect(configTextCommand("tree: t.nwk\n")).toStrictEqual(null);
  });
});

describe("config text with the current inputs", () => {
  const specs = commandSettings("timetree").specs;
  const current = setAt(setAt(defaultConfig(specs), ["tree"], "data/zika/86/tree.nwk"), ["alignment"], ["a.fasta"]);

  test("inputs the text does not name are appended", () => {
    expect(withMissingInputs("clock_rate: 0.0008\n", specs, current)).toStrictEqual(
      "clock_rate: 0.0008\nalignment:\n  - a.fasta\ntree: data/zika/86/tree.nwk\n",
    );
  });

  test("an input the text names is kept as written", () => {
    expect(withMissingInputs("tree: other.nwk\nalignment: []", specs, current)).toStrictEqual(
      "tree: other.nwk\nalignment: []",
    );
  });

  test("text that is not a mapping is left for the config check to report", () => {
    expect(withMissingInputs("- a\n- b\n", specs, current)).toStrictEqual("- a\n- b\n");
  });
});
