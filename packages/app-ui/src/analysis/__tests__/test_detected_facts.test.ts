import { commandSettings } from "../../settings/catalog";
import { detectedFacts } from "../detectedFacts";

const SPECS = commandSettings("timetree").settings;

const label = (key: string) => SPECS.find((spec) => spec.key === key)?.label;

describe("detected input facts of the data section", () => {
  test("each fact names the setting that changes it, and marks a value that a setting gives", () => {
    const facts = detectedFacts(
      {
        metadata_delimiter: { delimiter: "\t", origin: "detected", setting: "metadata_delimiters" },
        date_column: { column: "date", origin: "setting", setting: "date_column" },
        date_format: {
          readings: [
            { form: "calendar", format: "%Y-%m-%d", count: 98, text: "YYYY-MM-DD" },
            { form: "numeric", count: 2, text: "decimal years" },
          ],
          read: 100,
          rows: 100,
          origin: "detected",
          setting: "date_format",
        },
        sequence_length: { length: 18_959, origin: "detected", setting: "alignment" },
      },
      SPECS,
    );

    expect(facts).toStrictEqual([
      { label: label("metadata_delimiters"), value: "tab", setting: "metadata_delimiters", note: undefined },
      { label: label("date_column"), value: "column date", setting: "date_column", note: "set" },
      {
        label: label("date_format"),
        value: "100 of 100 read as YYYY-MM-DD, decimal years",
        setting: "date_format",
        note: undefined,
      },
      { label: "Sequence length", value: "18,959 sites", setting: "alignment", note: undefined },
    ]);
  });

  test("each fact takes the label of its setting from the catalog", () => {
    expect([label("metadata_delimiters"), label("date_column"), label("date_format")]).not.toContain(undefined);
  });

  test("no detected setting gives no fact", () => {
    expect(detectedFacts({}, SPECS)).toStrictEqual([]);
  });
});
