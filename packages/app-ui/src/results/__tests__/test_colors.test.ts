import { describe, expect, test } from "vitest";
import * as z from "zod";

import { mutedAuspiceDocument } from "../colors";
import { parseAuspiceJson, readAuspiceTree } from "../tree";

const zColorings = z.object({
  meta: z.object({ colorings: z.array(z.object({ key: z.string(), scale: z.unknown().optional() })) }),
});

describe("categorical colours", () => {
  test("a trait with few states takes Paul Tol's muted palette in alphabetical state order", () => {
    expect(scales(["south", "north", "south"])).toStrictEqual([
      [
        "region",
        [
          ["north", "#332288"],
          ["south", "#88ccee"],
        ],
      ],
      ["gt", undefined],
      ["num_date", undefined],
    ]);
  });

  test("a trait with more states than the nine palette colours keeps the Auspice default", () => {
    expect(scales(Array.from({ length: 10 }, (_, index) => `state${index}`))[0]).toStrictEqual(["region", undefined]);
  });

  test("nine states use the whole palette", () => {
    const scale = scales(Array.from({ length: 9 }, (_, index) => `state${index}`))[0]?.[1];

    expect(scale).toStrictEqual([
      ["state0", "#332288"],
      ["state1", "#88ccee"],
      ["state2", "#44aa99"],
      ["state3", "#117733"],
      ["state4", "#999933"],
      ["state5", "#ddcc77"],
      ["state6", "#cc6677"],
      ["state7", "#882255"],
      ["state8", "#aa4499"],
    ]);
  });
});

function scales(states: readonly string[]): Array<[string, unknown]> {
  const text = JSON.stringify({
    meta: {
      colorings: [
        { key: "region", title: "Region", type: "categorical" },
        { key: "gt", title: "Genotype", type: "categorical" },
        { key: "num_date", title: "Date", type: "continuous" },
      ],
    },
    tree: {
      name: "root",
      node_attrs: { region: { value: states[0] }, gt: { value: "A" } },
      children: states.map((state, index) => ({ name: `tip${index}`, node_attrs: { region: { value: state } } })),
    },
  });

  return zColorings
    .parse(mutedAuspiceDocument(text, readAuspiceTree(parseAuspiceJson(text))))
    .meta.colorings.map((coloring) => [coloring.key, coloring.scale]);
}
