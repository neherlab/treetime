import { detectedFacts } from "../detectedFacts";

describe("detected input facts of the data section", () => {
  test("each fact names the setting that changes it, and marks a value that a setting gives", () => {
    const facts = detectedFacts({
      metadata_delimiter: { delimiter: "\t", origin: "detected", setting: "metadata_delimiters" },
      date_column: { column: "date", origin: "setting", setting: "date_column" },
      date_format: {
        readings: [
          { form: "calendar", format: "%Y-%m-%d", count: 98 },
          { form: "numeric", count: 2 },
        ],
        read: 100,
        rows: 100,
        origin: "detected",
        setting: "date_format",
      },
      sequence_length: { length: 18_959, origin: "detected", setting: "alignment" },
    });

    expect(facts).toStrictEqual([
      { label: "Delimiter", value: "tab", setting: "metadata_delimiters", note: undefined },
      { label: "Dates", value: "column date", setting: "date_column", note: "set" },
      {
        label: "Date format",
        value: "100 of 100 read as %Y-%m-%d, decimal years",
        setting: "date_format",
        note: undefined,
      },
      { label: "Sequence length", value: "18,959 sites", setting: "alignment", note: undefined },
    ]);
  });

  test("no detected setting gives no fact", () => {
    expect(detectedFacts({})).toStrictEqual([]);
  });
});
