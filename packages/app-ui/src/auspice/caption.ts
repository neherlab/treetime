import escape from "escape-html";

import type { AuspicePublication } from "./state";

export function figureCaption({
  runs,
  version,
  updated,
  summary,
  mapAttribution,
  publications,
}: FigureCaptionParams): string[] {
  return [
    ...runLines(runs),
    version === undefined ? "Made with TreeTime" : `Made with TreeTime ${escape(version)}`,
    ...(updated === undefined ? [] : [`Updated ${escape(updated)}`]),
    escape(summary),
    ...(mapAttribution === undefined ? [] : [`Map: ${mapAttribution}`]),
    "",
    "Relevant publications:",
    ...publications.map(publicationLine),
  ];
}

export interface FigureCaptionParams {
  runs: readonly [CaptionRun] | readonly [CaptionRun, CaptionRun];
  version: string | undefined;
  updated: string | undefined;
  summary: string;
  mapAttribution: string | undefined;
  publications: readonly AuspicePublication[];
}

export interface CaptionRun {
  title: string;
  command: string;
}

function runLines(runs: FigureCaptionParams["runs"]): string[] {
  if (runs.length === 1) {
    const [run] = runs;

    return [escape(run.title), `Command: treetime ${escape(run.command)}`];
  }

  const [first, second] = runs;

  return [`First tree: ${runText(first)}`, `Second tree: ${runText(second)}`];
}

function publicationLine({ author, title, journal, year, href }: AuspicePublication): string {
  const text = `${author}, ${title}, ${journal} (${year})`;

  return `<a href="${escape(href)}">${escape(text)}</a>`;
}

function runText(run: CaptionRun): string {
  return `${escape(run.title)} (treetime ${escape(run.command)})`;
}
