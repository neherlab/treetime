import type { InputFacts, PreviewRoot, RegressionStatistic } from "@neherlab/app-contracts";

import { formatRate, yearDateText } from "../format";
import { baseName } from "../settings/inputs";

const STATISTIC_DIGITS = 3;

export const PREVIEW_ROOT_CAPTIONS: Readonly<Record<PreviewRoot, string>> = {
  regression: "Divergence from the root that the root-to-tip regression places",
  input: "Divergence from the root of the input tree, without a root-to-tip regression",
};

export function inputFactRows({ tree, alignments, metadata, tips_without_metadata }: InputFacts): FactRow[] {
  const rows: FactRow[] = [];

  if (tree !== undefined) {
    rows.push({ label: "Tree", value: treeText(tree.tips, tree.internal_nodes, tree.polytomies) });
  }

  for (const alignment of alignments) {
    rows.push({
      label: alignments.length === 1 ? "Alignment" : `Alignment ${alignment.number}`,
      value: `${baseName(alignment.path)}: ${count(alignment.sequences, "sequence")}, ${lengthText(alignment.min_length, alignment.max_length)}`,
    });
  }

  if (metadata !== undefined) {
    rows.push({ label: "Metadata", value: `${count(metadata.rows, "row")}, ID column ${metadata.id_column}` });

    if (metadata.dates !== undefined) {
      const { readable, unreadable, span } = metadata.dates;
      const range = span === undefined ? "" : `, ${yearDateText(span.earliest)} to ${yearDateText(span.latest)}`;
      const unread = unreadable.length === 0 ? "" : `, ${unreadable.length} unreadable`;

      rows.push({ label: "Dates", value: `${readable} readable${unread}${range}` });
    }
  }

  if (tips_without_metadata !== undefined && tips_without_metadata.length > 0) {
    rows.push({ label: "Tips without metadata", value: String(tips_without_metadata.length) });
  }

  return rows;
}

export function statisticText({ statistic, value, unit }: RegressionStatistic): string {
  if (value.kind === "date") {
    return yearDateText(value.date);
  }

  const number = statistic === "clock-rate" ? formatRate(value.value) : value.value.toPrecision(STATISTIC_DIGITS);

  return unit === undefined ? number : `${number} ${unit}`;
}

export interface FactRow {
  label: string;
  value: string;
}

function treeText(tips: number, internal: number, polytomies: number): string {
  const parts = [count(tips, "tip"), count(internal, "internal node")];

  return (polytomies > 0 ? [...parts, count(polytomies, "polytomy", "polytomies")] : parts).join(", ");
}

function lengthText(min: number, max: number): string {
  return min === max ? `length ${min}` : `length ${min} to ${max}`;
}

function count(value: number, singular: string, plural = `${singular}s`): string {
  return `${value} ${value === 1 ? singular : plural}`;
}
