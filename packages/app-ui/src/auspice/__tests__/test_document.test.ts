import { describe, expect, test } from "vitest";

import { withColorScales } from "../document";

describe("color scales in the Auspice document", () => {
  test("color scales reach the colorings they name and leave the original document unchanged", () => {
    const document = {
      meta: {
        colorings: [
          { key: "region", title: "Region", type: "categorical" },
          { key: "num_date", title: "Date", type: "continuous" },
        ],
      },
      tree: { name: "root" },
    };

    const scaled = withColorScales(
      document,
      new Map([
        [
          "region",
          [
            ["north", "#332288"],
            ["south", "#88ccee"],
          ],
        ],
      ]),
    );

    expect([scaled, document.meta.colorings[0]]).toStrictEqual([
      {
        meta: {
          colorings: [
            {
              key: "region",
              title: "Region",
              type: "categorical",
              scale: [
                ["north", "#332288"],
                ["south", "#88ccee"],
              ],
            },
            { key: "num_date", title: "Date", type: "continuous" },
          ],
        },
        tree: { name: "root" },
      },
      { key: "region", title: "Region", type: "categorical" },
    ]);
  });

  test("a document without colorings stays as it is", () => {
    const document = { meta: {}, tree: { name: "root" } };

    expect(withColorScales(document, new Map([["region", [["north", "#332288"]]]]))).toBe(document);
  });
});
