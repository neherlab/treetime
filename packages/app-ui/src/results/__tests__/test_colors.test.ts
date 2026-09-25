import { describe, expect, test } from "vitest";

import { mutedColorScales } from "../colors";
import type { ResultColoring } from "../types";

describe("categorical colours", () => {
  test("a trait with few states takes Paul Tol's muted palette in the order of its states", () => {
    expect([...mutedColorScales(colorings(["north", "south"]))]).toStrictEqual([
      [
        "region",
        [
          ["north", "#332288"],
          ["south", "#88ccee"],
        ],
      ],
    ]);
  });

  test("a trait with more states than the nine palette colours keeps the Auspice default", () => {
    expect(mutedColorScales(colorings(Array.from({ length: 10 }, (_, index) => `state${index}`))).size).toBe(0);
  });

  test("nine states use the whole palette", () => {
    const scale = mutedColorScales(colorings(Array.from({ length: 9 }, (_, index) => `state${index}`))).get("region");

    expect(scale?.map(([, color]) => color)).toStrictEqual([
      "#332288",
      "#88ccee",
      "#44aa99",
      "#117733",
      "#999933",
      "#ddcc77",
      "#cc6677",
      "#882255",
      "#aa4499",
    ]);
  });
});

function colorings(states: string[]): ResultColoring[] {
  return [
    { key: "region", title: "Region", kind: "categorical", states },
    { key: "gt", title: "Genotype", kind: "categorical", states: ["A"] },
    { key: "num_date", title: "Date", kind: "continuous", states: [] },
  ];
}
