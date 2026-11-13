import type {
  AlphabetName,
  DetectedDateFormat,
  DetectedSettings,
  SettingOrigin,
  SettingSpec,
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

const ORIGIN_NOTES: Readonly<Record<SettingOrigin, string | undefined>> = {
  detected: undefined,
  setting: "set",
  default: "default",
};

export function detectedFacts(detected: DetectedSettings, specs: readonly SettingSpec[]): DetectedFact[] {
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

  const fact = (
    { origin, setting }: { origin: SettingOrigin; setting: string },
    value: string,
    fileFact?: string,
  ): DetectedFact => {
    const spec = specs.find((candidate) => candidate.key === setting);
    const label = spec === undefined || spec.role === "input" ? (fileFact ?? spec?.label ?? setting) : spec.label;

    return { label, value, setting, note: ORIGIN_NOTES[origin] };
  };

  return [
    newick_dialect === undefined ? undefined : fact(newick_dialect, newick_dialect.dialect),
    metadata_delimiter === undefined
      ? undefined
      : fact(metadata_delimiter, DELIMITER_NAMES.get(metadata_delimiter.delimiter) ?? metadata_delimiter.delimiter),
    id_column === undefined ? undefined : fact(id_column, `column ${id_column.column}`),
    date_column === undefined ? undefined : fact(date_column, `column ${date_column.column}`),
    date_format === undefined ? undefined : fact(date_format, dateFormatText(date_format)),
    alphabet === undefined ? undefined : fact(alphabet, ALPHABET_NAMES[alphabet.alphabet]),
    sequence_length === undefined
      ? undefined
      : fact(sequence_length, `${sequence_length.length.toLocaleString("en-US")} sites`, "Sequence length"),
    partitions === undefined ? undefined : fact(partitions, partitionsText(partitions.partitions.length), "Partitions"),
  ].filter((entry) => entry !== undefined);
}

export interface DetectedFact {
  label: string;
  value: string;
  setting: string;
  note: string | undefined;
}

function dateFormatText({ readings, read, rows }: DetectedDateFormat): string {
  const forms = readings.map((reading) => reading.text).join(", ");

  return forms === "" ? `${read} of ${rows} read` : `${read} of ${rows} read as ${forms}`;
}

function partitionsText(count: number): string {
  return count === 1 ? "1 partition" : `${count} partitions`;
}
