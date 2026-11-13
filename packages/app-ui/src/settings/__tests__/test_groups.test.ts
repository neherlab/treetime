import { commandSettings } from "../catalog";
import { settingSections } from "../groups";

const SETTINGS = commandSettings("timetree");

describe("setting sections", () => {
  test("group the settings by purpose, the analysis settings by heading", () => {
    expect(
      settingSections(SETTINGS, SETTINGS.settings).map(({ title, groups }) => [
        title,
        groups.map(({ group }) => group),
      ]),
    ).toStrictEqual([
      ["Input data", ["Input data"]],
      [
        "Analysis",
        [
          "Molecular clock",
          "Rooting",
          "Branch lengths",
          "Dating",
          "Polytomies",
          "Coalescent prior",
          "Substitution model",
          "Ancestral reconstruction",
        ],
      ],
      ["Output format", ["Output"]],
      ["Tree ordering and reproducibility", ["Tree ordering", "Reproducibility"]],
    ]);
  });

  test("leave out the settings that the app sets for every run", () => {
    const output = settingSections(SETTINGS, SETTINGS.settings).find((section) => section.purpose === "output");

    expect(output?.groups.flatMap(({ specs }) => specs.map((spec) => spec.key))).toStrictEqual([
      "newick_annotations",
      "divergence_units",
      "ambiguity_threshold",
      "ambiguity_outputs",
    ]);
  });

  test("keep only the sections with a matching setting", () => {
    const seed = SETTINGS.settings.filter((spec) => spec.key === "seed");

    expect(settingSections(SETTINGS, seed).map(({ title }) => title)).toStrictEqual([
      "Tree ordering and reproducibility",
    ]);
  });
});
