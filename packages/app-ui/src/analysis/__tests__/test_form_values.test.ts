import { fromFormConfig, toFormConfig, toFormValue } from "../formValues";

describe("form values", () => {
  test("keeps scalars, lists and lists of records with scalar and list fields", () => {
    const config = {
      max_iter: 2,
      relax: [1, 0],
      entries: [{ name: "cp1", sites: ["101-1600\\3"], weight: 2, on: true }, {}],
    };

    expect(fromFormConfig(toFormConfig(config))).toStrictEqual(config);
  });

  test.each([
    ["a single object", { name: "cp1" }, 'The settings form cannot hold the value {"name":"cp1"}'],
    [
      "a record with an object field",
      [{ model: { name: "hky85" } }],
      'The settings form cannot hold the value [{"model":{"name":"hky85"}}]',
    ],
    [
      "a record with a list of lists",
      [{ sites: [["1-3"]] }],
      'The settings form cannot hold the value [{"sites":[["1-3"]]}]',
    ],
    [
      "records mixed with scalars",
      [{ name: "cp1" }, "cp2"],
      'The settings form cannot hold the value [{"name":"cp1"},"cp2"]',
    ],
  ])("rejects %s", (_, value, message) => {
    expect(() => toFormValue(value)).toThrow(new Error(message));
  });
});
