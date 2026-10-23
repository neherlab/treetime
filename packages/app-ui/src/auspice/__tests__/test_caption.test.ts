import { figureCaption } from "../caption";

const TREETIME = {
  author: "Sagulenko et al.",
  title: "TreeTime: Maximum-likelihood phylodynamic analysis",
  journal: "Virus Evolution",
  year: "2017",
  href: "https://doi.org/10.1093/ve/vex042",
};

const PUBLICATION_LINE =
  '<a href="https://doi.org/10.1093/ve/vex042">Sagulenko et al., TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017)</a>';

describe("figure caption", () => {
  test("one run names the run, the command, the version, the date and the publications", () => {
    expect(
      figureCaption({
        runs: [{ title: "Zika timetree", command: "timetree" }],
        version: "1.0.0",
        updated: "2026-10-24",
        summary: "Showing 86 of 86 genomes.",
        mapAttribution: undefined,
        publications: [TREETIME],
      }),
    ).toStrictEqual([
      "Zika timetree",
      "Command: treetime timetree",
      "Made with TreeTime 1.0.0",
      "Updated 2026-10-24",
      "Showing 86 of 86 genomes.",
      "",
      "Relevant publications:",
      PUBLICATION_LINE,
    ]);
  });

  test("two runs name the first and the second tree", () => {
    expect(
      figureCaption({
        runs: [
          { title: "Clock rate fixed", command: "timetree" },
          { title: "Clock rate inferred", command: "timetree" },
        ],
        version: "1.0.0",
        updated: undefined,
        summary: "Showing 20 of 20 genomes.",
        mapAttribution: undefined,
        publications: [],
      }),
    ).toStrictEqual([
      "First tree: Clock rate fixed (treetime timetree)",
      "Second tree: Clock rate inferred (treetime timetree)",
      "Made with TreeTime 1.0.0",
      "Showing 20 of 20 genomes.",
      "",
      "Relevant publications:",
    ]);
  });

  test("a figure with the map names the map data with its links", () => {
    const attribution = '<a href="https://www.naturalearthdata.com/">Natural Earth</a> 5.1.2';

    expect(
      figureCaption({
        runs: [{ title: "Countries", command: "mugration" }],
        version: undefined,
        updated: undefined,
        summary: "Showing 86 of 86 genomes.",
        mapAttribution: attribution,
        publications: [],
      }),
    ).toStrictEqual([
      "Countries",
      "Command: treetime mugration",
      "Made with TreeTime",
      "Showing 86 of 86 genomes.",
      `Map: ${attribution}`,
      "",
      "Relevant publications:",
    ]);
  });

  test("markup in a run title is written as text", () => {
    const [title, , , summary] = figureCaption({
      runs: [{ title: `<b>Flu</b> & "H3N2" 'all'`, command: "timetree" }],
      version: "1.0.0",
      updated: undefined,
      summary: "Showing <1 & 2>",
      mapAttribution: undefined,
      publications: [],
    });

    expect(title).toBe("&lt;b&gt;Flu&lt;/b&gt; &amp; &quot;H3N2&quot; &#39;all&#39;");
    expect(summary).toBe("Showing &lt;1 &amp; 2&gt;");
  });
});
