import { describe, expect, test } from "vitest";

import { commandLineLines, commandLineText, yamlLines, yamlText } from "../commandLine";
import { defaultConfig } from "../config";
import { setAt, type JsonObject, type JsonValue } from "../json";
import { commandSettings } from "../schema";

const timetree = commandSettings("timetree").specs;

const clock = commandSettings("clock").specs;

function withSettings(base: JsonObject, settings: ReadonlyArray<readonly [string[], JsonValue]>): JsonObject {
  let config = base;

  for (const [path, value] of settings) {
    config = setAt(config, path, value);
  }

  return config;
}

const RESOLVED_TIMETREE = withSettings(defaultConfig(timetree), [
  [["tree"], "data/zika/86/tree.nwk"],
  [["alignment"], ["data/zika/86/aln.fasta.xz"]],
  [["metadata"], "data/zika/86/metadata.tsv"],
  [["relax"], [1, 0]],
  [["model_params"], ["kappa=0.2", "pis=0.25,0.25,0.25,0.25"]],
  [["confidence"], true],
  [["reroot"], "least-squares"],
  [["output_selection"], ["Nwk", "Auspice", "Tracelog"]],
  [["output_all"], "out"],
]);

describe("command line", () => {
  test("writes inputs, changed settings in group order, the run additions and the output folder", () => {
    expect(commandLineText(commandLineLines("timetree", timetree, RESOLVED_TIMETREE))).toStrictEqual(
      [
        "treetime timetree",
        "  --tree data/zika/86/tree.nwk",
        "  --alignment data/zika/86/aln.fasta.xz",
        "  --metadata data/zika/86/metadata.tsv",
        "  --relax 1 0",
        "  --reroot least-squares",
        "  --confidence",
        "  --model-params kappa=0.2 --model-params pis=0.25,0.25,0.25,0.25",
        "  --output-selection nwk,auspice,tracelog",
        "  --output-all out",
      ].join(" \\\n"),
    );
  });

  test("marks the changed settings", () => {
    const kinds = commandLineLines("timetree", timetree, RESOLVED_TIMETREE).map((line) => line.kind);

    expect(kinds).toStrictEqual([
      "command",
      "input",
      "input",
      "input",
      "changed",
      "changed",
      "changed",
      "changed",
      "changed",
      "output",
    ]);
  });

  test("writes nested settings with their own flag", () => {
    const config = withSettings(defaultConfig(clock), [
      [["tree"], "t.nwk"],
      [["branch_split", "n_points"], 21],
      [["clock_regression", "variance_factor"], 0.5],
    ]);

    expect(commandLineLines("clock", clock, config).map((line) => line.text)).toStrictEqual([
      "treetime clock",
      "--tree t.nwk",
      "--variance-factor 0.5",
      "--branch-split-grid-n-points 21",
      "--output-all out",
    ]);
  });

  test("quotes values the shell would split", () => {
    const config = withSettings(defaultConfig(clock), [[["tree"], "my data/tree one.nwk"]]);

    expect(commandLineLines("clock", clock, config)[1]?.text).toStrictEqual("--tree 'my data/tree one.nwk'");
  });

  test("a value without a command-line form becomes a comment", () => {
    const config = withSettings(defaultConfig(timetree), [[["metadata_id_columns"], []]]);

    expect(commandLineLines("timetree", timetree, config)[1]).toStrictEqual({
      text: "# metadata_id_columns = [] has no command-line form; use the YAML config",
      kind: "comment",
    });
  });

  test("the command-line text leaves comments out", () => {
    const config = withSettings(defaultConfig(timetree), [[["metadata_id_columns"], []]]);

    expect(commandLineText(commandLineLines("timetree", timetree, config))).toStrictEqual(
      "treetime timetree \\\n  --output-all out",
    );
  });
});

describe("yaml config", () => {
  test("writes inputs, changed settings, the run additions and the output folder", () => {
    expect(yamlText(yamlLines("timetree", timetree, RESOLVED_TIMETREE))).toStrictEqual(
      [
        "# yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-timetree.schema.json",
        "# treetime timetree --config run.yaml",
        "tree: data/zika/86/tree.nwk",
        "alignment:",
        "  - data/zika/86/aln.fasta.xz",
        "metadata: data/zika/86/metadata.tsv",
        "relax:",
        "  - 1",
        "  - 0",
        "reroot: least-squares",
        "confidence: true",
        "model_params:",
        "  - kappa=0.2",
        "  - pis=0.25,0.25,0.25,0.25",
        "output_selection:",
        "  - Nwk",
        "  - Auspice",
        "  - Tracelog",
        "output_all: out",
        "",
      ].join("\n"),
    );
  });

  test("writes nested settings under their parent key", () => {
    const config = withSettings(defaultConfig(clock), [
      [["branch_split", "n_points"], 21],
      [["branch_split", "method"], "brent"],
    ]);

    expect(yamlLines("clock", clock, config).slice(2)).toStrictEqual([
      { text: "branch_split:", kind: "changed" },
      { text: "  method: brent", kind: "changed" },
      { text: "  n_points: 21", kind: "changed" },
      { text: "output_all: out", kind: "output" },
    ]);
  });
});
