import type {
  DrmAnnotation,
  GapFill,
  HomoplasySite,
  HomoplasyStatistics,
  MultiplicityRow,
  RecurrentMutation,
  RunRecord,
  SiteHitsRow,
} from "@neherlab/app-contracts";

import { nucleotideColorBy } from "../auspice/genotype";
import { runTimeEntry, type SummaryEntry } from "../components/Panel";
import type { AxisFrame } from "./palette";

const LOG_AXIS_FLOOR = 0.01;

const GAP_FILL_NOTES: Readonly<Record<GapFill, string | undefined>> = {
  "only-terminal":
    "Leading and trailing gaps became N before reconstruction, so these counts include missing coverage at the sequence ends.",
  all: "Gaps became N before reconstruction, so these counts include missing coverage.",
  none: undefined,
};

export interface SiteHitsPoint extends SiteHitsRow {
  sitesShown: number | undefined;
  expectedShown: number | undefined;
}

export interface MultiplicityPoint extends MultiplicityRow {
  mutationsShown: number | undefined;
}

export function initialHomoplasyColorBy(statistics: HomoplasyStatistics | undefined): string | undefined {
  const top = statistics?.recurrent[0];

  return top === undefined ? undefined : nucleotideColorBy(top.position);
}

export function siteAt(sites: readonly HomoplasySite[], position: number | undefined): HomoplasySite | undefined {
  return position === undefined ? undefined : sites.find((site) => site.position === position);
}

export function pressedMutations(rows: readonly RecurrentMutation[], position: number | undefined): Set<string> {
  return new Set(rows.flatMap((row) => (row.position === position ? [row.mutation] : [])));
}

export function homoplasySummary(record: RunRecord, statistics: HomoplasyStatistics): SummaryEntry[] {
  const drm = statistics.recurrent_drm_substitutions;
  const drmText = drm === undefined ? "" : `, ${drm} at drug resistance positions`;

  return [
    {
      label: "Recurrent substitutions",
      value: String(statistics.recurrent_substitutions),
      detail: `of ${statistics.distinct_substitutions} distinct substitutions${drmText}`,
    },
    {
      label: "Sites hit more than once",
      value: String(statistics.sites_hit_more_than_once),
      detail: `Poisson expectation ${statistics.expected_sites_hit_more_than_once.toFixed(1)}`,
    },
    {
      label: "Poisson log-likelihood difference",
      value: statistics.log_likelihood_difference.toExponential(2),
      detail: "Negative: substitutions cluster at fewer sites than expected",
    },
    {
      label: "Samples with homoplasies",
      value: String(statistics.samples_with_homoplasies),
      detail: "Terminal branch has a substitution at a site hit more than once",
    },
    {
      label: "Ambiguous changes",
      value: String(statistics.ambiguous_changes),
      detail: "Changes to or from N and other ambiguity codes",
    },
    runTimeEntry(record),
  ];
}

export function gapFillNote(gapFill: GapFill | undefined): string | undefined {
  return gapFill === undefined ? undefined : GAP_FILL_NOTES[gapFill];
}

export function drmText(drm: DrmAnnotation): string {
  return [drm.gene, drm.drug, drm.substitution].filter((part) => part !== undefined).join(" ");
}

export function siteText(site: HomoplasySite): string {
  return `Position ${site.display_position}, ${site.branches} branches`;
}

export function logAxis(values: readonly number[]): AxisFrame {
  const top = Math.ceil(Math.log10(Math.max(1, ...values)));
  const bottom = Math.log10(LOG_AXIS_FLOOR);
  const ticks = Array.from({ length: top - bottom + 1 }, (_, index) => 10 ** (bottom + index));

  return { domain: [LOG_AXIS_FLOOR, 10 ** top], ticks };
}

export function siteHitsPoints(rows: readonly SiteHitsRow[]): SiteHitsPoint[] {
  return rows.map((row) => ({ ...row, sitesShown: onLogAxis(row.sites), expectedShown: onLogAxis(row.expected) }));
}

export function multiplicityPoints(rows: readonly MultiplicityRow[]): MultiplicityPoint[] {
  return rows.map((row) => ({ ...row, mutationsShown: onLogAxis(row.mutations) }));
}

function onLogAxis(value: number): number | undefined {
  return value >= LOG_AXIS_FLOOR ? value : undefined;
}
