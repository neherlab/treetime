import type {
  AlphabetName,
  DateForm,
  DateReading,
  DetectedDateFormat,
  DetectedSettings,
  SettingOrigin,
} from "@neherlab/app-contracts";

const DELIMITER_NAMES = new Map([
  ["\t", "tab"],
  [",", "comma"],
  [";", "semicolon"],
  [" ", "space"],
  ["|", "vertical bar"],
]);

const ALPHABET_NAMES: Readonly<Record<AlphabetName, string>> = {
  nuc: "nucleotides",
  aa: "amino acids",
  "aa-no-stop": "amino acids without stop codon",
};

const READING_TEXT: Readonly<Record<DateForm, string>> = {
  numeric: "decimal years",
  calendar: "calendar dates",
  uncertain: "uncertain dates",
  range: "date ranges",
};

const ORIGIN_NOTES: Readonly<Record<SettingOrigin, string | undefined>> = {
  detected: undefined,
  setting: "set",
  default: "default",
};

export function detectedFacts(detected: DetectedSettings): DetectedFact[] {
  const {
    newick_dialect,
    metadata_delimiter,
    id_column,
    date_column,
    date_format,
    alphabet,
    sequence_length,
    partitions,
  } = detected;

  return [
    newick_dialect === undefined ? undefined : fact("Newick dialect", newick_dialect, newick_dialect.dialect),
    metadata_delimiter === undefined
      ? undefined
      : fact(
          "Delimiter",
          metadata_delimiter,
          DELIMITER_NAMES.get(metadata_delimiter.delimiter) ?? metadata_delimiter.delimiter,
        ),
    id_column === undefined ? undefined : fact("Sample names", id_column, `column ${id_column.column}`),
    date_column === undefined ? undefined : fact("Dates", date_column, `column ${date_column.column}`),
    date_format === undefined ? undefined : fact("Date format", date_format, dateFormatText(date_format)),
    alphabet === undefined ? undefined : fact("Alphabet", alphabet, ALPHABET_NAMES[alphabet.alphabet]),
    sequence_length === undefined
      ? undefined
      : fact("Sequence length", sequence_length, `${sequence_length.length.toLocaleString("en-US")} sites`),
    partitions === undefined ? undefined : fact("Partitions", partitions, partitionsText(partitions.partitions.length)),
  ].filter((entry) => entry !== undefined);
}

export interface DetectedFact {
  label: string;
  value: string;
  setting: string;
  note: string | undefined;
}

function fact(
  label: string,
  { origin, setting }: { origin: SettingOrigin; setting: string },
  value: string,
): DetectedFact {
  return { label, value, setting, note: ORIGIN_NOTES[origin] };
}

function dateFormatText({ readings, read, rows }: DetectedDateFormat): string {
  const forms = readings.map(readingText).join(", ");

  return forms === "" ? `${read} of ${rows} read` : `${read} of ${rows} read as ${forms}`;
}

function readingText(reading: DateReading): string {
  return reading.form === "calendar" ? (reading.format ?? READING_TEXT.calendar) : READING_TEXT[reading.form];
}

function partitionsText(count: number): string {
  return count === 1 ? "1 partition" : `${count} partitions`;
}
