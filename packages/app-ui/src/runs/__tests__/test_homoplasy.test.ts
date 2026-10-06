import type { HomoplasyStatistics, RecurrentMutation, RunRecord } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";

import {
  drmText,
  gapFillNote,
  homoplasySummary,
  initialHomoplasyColorBy,
  logAxis,
  pressedMutations,
  siteAt,
  siteHitsPoints,
  siteText,
} from "../homoplasy";

const G10T: RecurrentMutation = {
  mutation: "G10T",
  position: 10,
  display_position: 10,
  branches: 2,
  terminal_branches: 2,
  branch_names: ["C", "D"],
};

const G10A: RecurrentMutation = {
  ...G10T,
  mutation: "G10A",
  branches: 3,
  branch_names: ["A", "B", "node_1"],
  drm: { gene: "HA", drug: "oseltamivir", substitution: "H275Y" },
};

const T20C: RecurrentMutation = { ...G10T, mutation: "T20C", position: 20, display_position: 20 };

const WITHOUT_DRMS: HomoplasyStatistics = {
  drm_annotated: false,
  zero_based: false,
  genome_length: 100,
  total_branch_length: 0.5,
  substitutions: 8,
  distinct_substitutions: 4,
  recurrent_substitutions: 3,
  sites_hit_more_than_once: 2,
  expected_sites_hit_more_than_once: 0.234,
  log_likelihood_difference: -3.5,
  samples_with_homoplasies: 1,
  ambiguous_changes: 4,
  indels: 3,
  site_hits: [],
  multiplicities: [],
  recurrent: [G10A, G10T, T20C],
  sites: [
    { position: 10, display_position: 10, branches: 5, substitutions: [] },
    { position: 20, display_position: 20, branches: 2, substitutions: [] },
  ],
  recurrent_indels: [],
  taxa: [],
  ambiguous_sites: [],
  ambiguous_site_count: 0,
};

const STATISTICS: HomoplasyStatistics = { ...WITHOUT_DRMS, drm_annotated: true, recurrent_drm_substitutions: 1 };

const RECORD: RunRecord = {
  id: "r1",
  title: "Homoplasy",
  command: "homoplasy",
  config: { tree: "t.nwk" },
  status: "ok",
  pinned: false,
  created_at: "2026-09-25T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [],
  changed_settings: [],
  headline: {},
  output_files: [],
  duration_seconds: 2.5,
  warnings: [],
};

describe("homoplasy results", () => {
  test("the tree opens colored by the position of the top recurrent mutation", () => {
    expect([
      initialHomoplasyColorBy(STATISTICS),
      initialHomoplasyColorBy({ ...STATISTICS, recurrent: [] }),
      initialHomoplasyColorBy(undefined),
    ]).toStrictEqual(["gt-nuc_10", undefined, undefined]);
  });

  test.each([
    { name: "site hit more than once", position: 20, expected: 20 },
    { name: "position without a site row", position: 30, expected: undefined },
    { name: "no genotype position", position: undefined, expected: undefined },
  ])("site at $name", ({ position, expected }) => {
    expect(siteAt(STATISTICS.sites, position)?.position).toStrictEqual(expected);
  });

  test.each([
    { name: "two alleles at one position", position: 10, expected: ["G10A", "G10T"] },
    { name: "one allele", position: 20, expected: ["T20C"] },
    { name: "no genotype position", position: undefined, expected: [] },
  ])("pressed mutations: $name", ({ position, expected }) => {
    expect([...pressedMutations(STATISTICS.recurrent, position)]).toStrictEqual(expected);
  });

  test("the summary counts recurrent substitutions, sites, samples, and ambiguous changes", () => {
    expect(homoplasySummary(RECORD, STATISTICS)).toStrictEqual([
      {
        label: "Recurrent substitutions",
        value: "3",
        detail: "of 4 distinct substitutions, 1 at drug resistance positions",
      },
      { label: "Sites hit more than once", value: "2", detail: "Poisson expectation 0.2" },
      {
        label: "Poisson log-likelihood difference",
        value: "-3.50e+0",
        detail: "Negative: substitutions cluster at fewer sites than expected",
      },
      {
        label: "Samples with homoplasies",
        value: "1",
        detail: "Terminal branch has a substitution at a site hit more than once",
      },
      { label: "Ambiguous changes", value: "4", detail: "Changes to or from N and other ambiguity codes" },
      { label: "Run time", value: "2.5 s" },
    ]);
  });

  test("the summary leaves out drug resistance without a DRM table", () => {
    const summary = homoplasySummary(RECORD, WITHOUT_DRMS);

    expect(summary[0]?.detail).toStrictEqual("of 4 distinct substitutions");
  });

  test.each([
    { gapFill: "only-terminal", expected: "at the sequence ends." },
    { gapFill: "all", expected: "include missing coverage." },
  ] as const)("gap filling $gapFill explains the ambiguous counts", ({ gapFill, expected }) => {
    expect(gapFillNote(gapFill)?.endsWith(expected)).toStrictEqual(true);
  });

  test("no gap filling needs no explanation", () => {
    expect([gapFillNote("none"), gapFillNote(undefined)]).toStrictEqual([undefined, undefined]);
  });

  test("a DRM annotation lists gene, drug, and substitution when known", () => {
    expect([
      drmText({ gene: "HA", drug: "oseltamivir", substitution: "H275Y" }),
      drmText({ gene: "RT", drug: "NRTI" }),
    ]).toStrictEqual(["HA oseltamivir H275Y", "RT NRTI"]);
  });

  test("a site reads as its display position and branch count", () => {
    expect(siteText({ position: 11, display_position: 10, branches: 4, substitutions: [] })).toStrictEqual(
      "Position 10, 4 branches",
    );
  });

  test("values below the log axis floor, zero included, are not drawn", () => {
    expect(
      siteHitsPoints([
        { hits: 0, sites: 90, expected: 93.2 },
        { hits: 5, sites: 0, expected: 0.004 },
      ]).map((point) => [point.sitesShown, point.expectedShown]),
    ).toStrictEqual([
      [90, 93.2],
      [undefined, undefined],
    ]);
  });

  test("the log axis spans powers of ten from the floor to above the largest value", () => {
    expect(logAxis([29_000, 0.5])).toStrictEqual({
      domain: [0.01, 100_000],
      ticks: [0.01, 0.1, 1, 10, 100, 1000, 10_000, 100_000],
    });
  });
});
