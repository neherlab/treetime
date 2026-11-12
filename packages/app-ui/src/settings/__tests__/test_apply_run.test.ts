import { withRunSettings } from "../applyRun";
import { commandSettings } from "../catalog";
import { normalizeConfig } from "../config";

const SPECS = commandSettings("timetree").settings;

describe("settings of a run applied to a draft", () => {
  test("only the settings that the run changed apply; the other settings and the inputs of the draft stay", () => {
    const draft = { tree: "draft/tree.nwk", clock_rate: 0.001, max_iter: 4 };

    const run = {
      tree: "run/tree.nwk",
      alignment: ["run/aln.fasta"],
      clock_rate: 0.002,
      coalescent: "skyline",
      confidence: true,
      max_iter: 2,
    };

    expect(withRunSettings(SPECS, draft, ["tree", "alignment", "coalescent", "confidence"], run)).toStrictEqual(
      expect.objectContaining({
        tree: "draft/tree.nwk",
        alignment: [],
        clock_rate: 0.001,
        max_iter: 4,
        coalescent: "skyline",
        confidence: true,
      }),
    );
  });

  test("the value of the run replaces the value of the draft", () => {
    const applied = withRunSettings(SPECS, { tree: "t.nwk", clock_filter: 3 }, ["clock_filter"], { clock_filter: 5 });

    expect([applied["tree"], applied["clock_filter"]]).toStrictEqual(["t.nwk", 5]);
  });

  test("settings that the app sets never apply", () => {
    const applied = withRunSettings(SPECS, { tree: "t.nwk" }, ["output_selection", "output_all"], {
      output_selection: ["nwk"],
      output_all: "out",
    });

    expect([applied["output_selection"], applied["output_all"]]).toStrictEqual([
      normalizeConfig(SPECS, { tree: "t.nwk" })["output_selection"],
      undefined,
    ]);
  });
});
