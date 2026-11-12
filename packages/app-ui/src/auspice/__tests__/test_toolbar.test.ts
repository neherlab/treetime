import { colorOptions, genotypeColorBy, layoutOptions, tipNames } from "../toolbar";

const COLORINGS = {
  num_date: { title: "Date", type: "continuous" },
  bad_branch: { title: "Excluded", type: "categorical" },
  gt: { title: "Genotype", type: "categorical" },
};

describe("viewer toolbar", () => {
  test("the color choices list the colorings of the tree without the bare genotype entry", () => {
    expect(colorOptions(COLORINGS, "num_date")).toStrictEqual([
      { value: "num_date", label: "Date" },
      { value: "bad_branch", label: "Excluded" },
    ]);
  });

  test("a genotype coloring in use is the first choice", () => {
    expect(colorOptions(COLORINGS, genotypeColorBy(123)).at(0)).toStrictEqual({
      value: "gt-nuc_123",
      label: "Genotype nuc 123",
    });
  });

  test("the clock layout needs both divergence and dates", () => {
    expect([layoutOptions(["div"]).length, layoutOptions(["div", "num_date"]).at(-1)]).toStrictEqual([
      3,
      { value: "clock", label: "Clock" },
    ]);
  });

  test("the search offers the sorted names of the tips", () => {
    expect(
      tipNames([
        { name: "NODE_1", hasChildren: true },
        { name: "B", hasChildren: false },
        { name: "A", hasChildren: false },
      ]),
    ).toStrictEqual(["A", "B"]);
  });
});
