import { describe, expect, test } from "vitest";

import { nucleotideColorBy, nucleotidePosition } from "../genotype";

describe("genotype color-by keys", () => {
  test("a nucleotide position encodes as Auspice's gt-nuc key", () => {
    expect(nucleotideColorBy(13_427)).toStrictEqual("gt-nuc_13427");
  });

  test.each([
    { name: "nucleotide position", colorBy: "gt-nuc_13427", expected: 13_427 },
    { name: "round trip", colorBy: nucleotideColorBy(5) ?? "", expected: 5 },
    { name: "metadata coloring", colorBy: "country", expected: undefined },
    { name: "amino-acid position", colorBy: "gt-HA1_144", expected: undefined },
    { name: "several positions", colorBy: "gt-nuc_142,144", expected: undefined },
  ])("$name", ({ colorBy, expected }) => {
    expect(nucleotidePosition(colorBy)).toStrictEqual(expected);
  });
});
