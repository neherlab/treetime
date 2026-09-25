import { readdirSync, readFileSync, statSync } from "node:fs";

import type { AppCommand } from "@neherlab/app-contracts";
import { tsvParse } from "d3-dsv";
import { describe, expect, test } from "vitest";
import * as z from "zod";

import { indexClades, matchAncestors } from "../clades";
import { mutedAuspiceDocument } from "../colors";
import { ancestorShifts, compareEstimates, timetreeEstimates } from "../estimates";
import type { RunFileEntry } from "../files";
import { loadRunResults } from "../load";
import { recurrentSites, stateChanges } from "../mutations";
import {
  readAugurClock,
  readClockCsv,
  readClockModel,
  readCoalescentTsv,
  readGtr,
  readTotalBranchLength,
  readTraitsCsv,
  readTracelog,
} from "../readers";
import { parseAuspiceJson, readAuspiceTree, type ResultTree } from "../tree";

const FIXTURE_ROOT = new URL("../../../../../tmp/app-fixtures/", import.meta.url);

const DATA_ROOT = new URL("../../../../../data/zika/86/", import.meta.url);

const NEWICK_TIP = /[(,]'(?<name>[^']+)':/gu;

const DATASET_TIPS = [...datasetText("tree.nwk").matchAll(NEWICK_TIP)].map((match) => match.groups?.["name"] ?? "");

const METADATA = tsvParse(datasetText("metadata.tsv"));

const KINDS: ReadonlyArray<[RegExp, RunFileEntry["kind"]]> = [
  [/\.auspice\.json$/u, "auspice"],
  [/\.clock-model\.json$/u, "clock-model"],
  [/\.augur-node-data\.json$/u, "augur-node-data"],
  [/\.gtr\.json$/u, "gtr"],
  [/\.tracelog\.csv$/u, "tracelog"],
  [/\.coalescent\.tsv$/u, "coalescent-tsv"],
  [/\.clock\.csv$/u, "clock-csv"],
  [/\.traits\.csv$/u, "traits-csv"],
];

describe("the tree reader on the zika/86 outputs", () => {
  test("the timetree tree has the samples of the input tree", () => {
    const tree = fixtureTree("timetree/timetree.auspice.json");

    expect(tree.tips.map((tip) => tip.name).toSorted(byText)).toStrictEqual(DATASET_TIPS.toSorted(byText));
  });

  test("every timetree node date lies inside its 90% interval", () => {
    const tree = fixtureTree("timetree/timetree.auspice.json");

    const outside = tree.nodes.filter(
      (node) =>
        node.date === undefined ||
        node.dateInterval === undefined ||
        node.date < node.dateInterval[0] ||
        node.date > node.dateInterval[1],
    );

    expect(outside.map((node) => node.name)).toStrictEqual([]);
  });

  test("node dates agree with the confidence table the same run wrote", () => {
    const tree = fixtureTree("timetree/timetree.auspice.json");

    const table = new Map(
      tsvParse(fixtureText("timetree/timetree.confidence.tsv")).map((row) => [row["name"], Number(row["date"])]),
    );

    const differences = tree.nodes
      .filter((node) => node !== tree.root)
      .map((node) => Math.abs((node.date ?? Number.NaN) - (table.get(node.name) ?? Number.NaN)));

    expect(differences.length).toBe(tree.nodes.length - 1);
    expect(Math.max(...differences)).toBeLessThanOrEqual(5e-4);
  });

  test("the mugration tree carries the trait of every sample from the metadata", () => {
    const tree = fixtureTree("mugration/mugration.auspice.json");
    const metadata = new Map(METADATA.map((row) => [row["#name"], row["country"]]));
    const mismatches = tree.tips.filter((tip) => tip.traits.get("country")?.value !== metadata.get(tip.name));

    expect(mismatches.map((tip) => tip.name)).toStrictEqual([]);
  });

  test("branch mutations match the node data of the same ancestral run", () => {
    const tree = fixtureTree("ancestral/ancestral.auspice.json");
    const nodeData = fixtureJson("ancestral/ancestral.augur-node-data.json");
    const expected = nodeMutations(nodeData);

    expect(new Map(tree.nodes.map((node) => [node.name, node.mutations.toSorted(byText)]))).toStrictEqual(expected);
  });

  test("prune and optimize trees read without dates or traits", () => {
    const trees = ["prune/prune.auspice.json", "optimize/optimize.auspice.json"].map((path) => fixtureTree(path));

    expect(
      trees.map((tree) => tree.nodes.some((node) => node.date !== undefined || node.traits.size > 0)),
    ).toStrictEqual([false, false]);
  });
});

describe("the clock readers on the zika/86 outputs", () => {
  test("the clock model rate is the rate of the node data", () => {
    const model = readClockModel(fixtureText("timetree/timetree.clock-model.json"));
    const augur = readAugurClock(fixtureText("timetree/timetree.augur-node-data.json"));

    expect(model.rate).toBe(augur?.rate);
  });

  test("an estimated clock has a correlation coefficient within [-1, 1]", () => {
    const model = readClockModel(fixtureText("timetree/timetree.clock-model.json"));

    expect({ fixed: model.fixed, bounded: model.r !== undefined && Math.abs(model.r) <= 1 }).toStrictEqual({
      fixed: false,
      bounded: true,
    });
  });

  test("the rate standard deviation is the square root of the rate variance", () => {
    const augur = readAugurClock(fixtureText("timetree/timetree.augur-node-data.json"));
    const variance = rateVariance(fixtureJson("timetree/timetree.clock-model.json"));

    expect(augur?.rateStd).toBeCloseTo(Math.sqrt(variance), 12);
  });

  test("the clock table has one row per node and flags outliers only on dated samples", () => {
    const rows = readClockCsv(fixtureText("clock/clock.clock.csv"));
    const tree = fixtureTree("clock/clock.auspice.json");
    const tips = new Set(tree.tips.map((tip) => tip.name));
    const outliers = rows.filter((row) => row.outlier);

    expect({
      rows: rows.length,
      someOutliers: outliers.length > 0,
      outliersAreDatedTips: outliers.every((row) => tips.has(row.name) && row.date !== undefined),
    }).toStrictEqual({ rows: tree.nodes.length, someOutliers: true, outliersAreDatedTips: true });
  });

  test("clock outliers are excluded tips in the Auspice tree of the same run", () => {
    const rows = readClockCsv(fixtureText("clock/clock.clock.csv"));
    const tree = fixtureTree("clock/clock.auspice.json");
    const excluded = new Set(tree.tips.filter((tip) => tip.excluded === true).map((tip) => tip.name));

    expect(rows.filter((row) => row.outlier && !excluded.has(row.name))).toStrictEqual([]);
  });
});

describe("the iteration and coalescent readers on the zika/86 outputs", () => {
  test("the tracelog has one row per iteration line", () => {
    const text = fixtureText("timetree/timetree.tracelog.csv");
    const lines = text.trim().split("\n").length - 1;

    expect(readTracelog(text).map((row) => row.iteration)).toStrictEqual([...Array.from({ length: lines }).keys()]);
  });

  test("an empty coalescent column reads as absent without a coalescent prior", () => {
    const rows = readTracelog(fixtureText("timetree/timetree.tracelog.csv"));

    expect(rows.map((row) => row.logLhCoal)).toStrictEqual(rows.map(() => undefined));
  });

  test("the total log likelihood is the sum of its parts under a skyline prior", () => {
    const rows = readTracelog(fixtureText("timetree-skyline/timetree.tracelog.csv"));

    const residuals = rows.map((row) =>
      Math.abs((row.logLhSeq ?? 0) + (row.logLhPos ?? 0) + (row.logLhCoal ?? 0) - (row.logLhTotal ?? Number.NaN)),
    );

    expect(Math.max(...residuals)).toBeLessThan(1e-6);
  });

  test("skyline segments tile the time axis and each estimate lies inside its interval", () => {
    const segments = readCoalescentTsv(fixtureText("timetree-skyline/timetree.coalescent.tsv"));
    const gaps = segments.slice(1).map((segment, index) => segment.start - (segments[index]?.end ?? Number.NaN));

    const outside = segments.filter(
      (segment) =>
        segment.ne.value < segment.ne.lower ||
        segment.ne.value > segment.ne.upper ||
        segment.tc.value < segment.tc.lower ||
        segment.tc.value > segment.tc.upper,
    );

    expect({ segments: segments.length, maxGap: Math.max(...gaps.map(Math.abs)), outside }).toStrictEqual({
      segments: 20,
      maxGap: 0,
      outside: [],
    });
  });

  test("the skyline table agrees with the coalescent JSON of the same run", () => {
    const segments = readCoalescentTsv(fixtureText("timetree-skyline/timetree.coalescent.tsv"));
    const json = skylineJsonValues(fixtureJson("timetree-skyline/timetree.coalescent.json"));

    expect(segments.map((segment) => [segment.start, segment.end, segment.ne.value])).toStrictEqual(json);
  });
});

describe("the trait and optimize readers on the zika/86 outputs", () => {
  test("the traits table names the attribute and agrees with the Auspice tree", () => {
    const traits = readTraitsCsv(fixtureText("mugration/mugration.traits.csv"));
    const tree = fixtureTree("mugration/mugration.auspice.json");
    const disagreeing = tree.nodes.filter((node) => node.traits.get("country")?.value !== traits.states.get(node.name));

    expect({ attribute: traits.attribute, disagreeing: disagreeing.map((node) => node.name) }).toStrictEqual({
      attribute: "country",
      disagreeing: [],
    });
  });

  test("the total branch length is the sum of the node data branch lengths", () => {
    const text = fixtureText("optimize/optimize.augur-node-data.json");

    const { nodes } = z
      .object({ nodes: z.record(z.string(), z.object({ branch_length: z.number() })) })
      .parse(JSON.parse(text));

    const expected = Object.values(nodes).reduce((sum, node) => sum + node.branch_length, 0);

    expect(readTotalBranchLength(text)).toBeCloseTo(expected, 12);
  });

  test("the substitution model names the model the run inferred", () => {
    expect(readGtr(fixtureText("optimize/optimize.gtr.json")).model).toBe("infer");
  });

  test("state changes in a real run agree with the traits table", () => {
    const tree = fixtureTree("mugration/mugration.auspice.json");
    const traits = readTraitsCsv(fixtureText("mugration/mugration.traits.csv"));

    const expected = [...tree.parents].filter(
      ([child, parent]) => traits.states.get(child.name) !== traits.states.get(parent.name),
    ).length;

    expect(stateChanges(tree, "country").reduce((sum, change) => sum + change.branches, 0)).toBe(expected);
  });

  test("the recurrent sites of a real run count each branch once", () => {
    const tree = fixtureTree("ancestral/ancestral.auspice.json");
    const total = recurrentSites(tree).reduce((sum, site) => sum + site.branches, 0);
    const bound = tree.nodes.reduce((sum, node) => sum + node.mutations.length, 0);

    expect(total <= bound && recurrentSites(tree).every((site) => site.branches > 1)).toBe(true);
  });
});

describe("comparison of the zika/86 time trees", () => {
  test("the roots of two runs on the same samples match", () => {
    const tree = fixtureTree("timetree/timetree.auspice.json");
    const skyline = fixtureTree("timetree-skyline/timetree.auspice.json");
    const matched = matchAncestors(indexClades(tree), indexClades(skyline));

    expect(matched.some((pair) => pair.first === tree.root && pair.second === skyline.root)).toBe(true);
  });

  test("the root shift between two real runs equals the shift of their matched roots", () => {
    const outputs = ["timetree", "timetree-skyline"].map((run) => ({
      tree: fixtureTree(`${run}/timetree.auspice.json`),
      clockModel: readClockModel(fixtureText(`${run}/timetree.clock-model.json`)),
      augurClock: readAugurClock(fixtureText(`${run}/timetree.augur-node-data.json`)),
      trace: readTracelog(fixtureText(`${run}/timetree.tracelog.csv`)),
    }));

    const [first, second] = outputs.map(timetreeEstimates);
    const [firstTree, secondTree] = outputs.map((output) => output.tree);

    if (first === undefined || second === undefined || firstTree === undefined || secondTree === undefined) {
      throw new Error("two runs expected");
    }

    const rootShift = ancestorShifts(firstTree, secondTree).find((shift) => shift.name === firstTree.root.name);

    expect(compareEstimates(first, second).rootShiftDays).toBe(rootShift?.shiftDays);
  });

  test("the summary of the zika run counts every sample and reads the final log likelihood", () => {
    const tree = fixtureTree("timetree/timetree.auspice.json");
    const trace = readTracelog(fixtureText("timetree/timetree.tracelog.csv"));

    const summary = timetreeEstimates({
      tree,
      clockModel: readClockModel(fixtureText("timetree/timetree.clock-model.json")),
      augurClock: readAugurClock(fixtureText("timetree/timetree.augur-node-data.json")),
      trace,
    });

    expect({
      samples: summary.samples,
      rootDate: summary.rootDate,
      rootInterval: summary.rootInterval,
      logLikelihood: summary.logLikelihood,
      iterations: summary.iterations,
    }).toStrictEqual({
      samples: 86,
      rootDate: tree.root.date,
      rootInterval: tree.root.dateInterval,
      logLikelihood: trace.at(-1)?.logLhTotal,
      iterations: trace.length,
    });
  });
});

describe("loading run results from the zika/86 outputs", () => {
  test("a skyline time tree run loads every output its view shows", async () => {
    const results = await load("timetree", "timetree-skyline");

    expect({
      tree: results.auspice?.tree.tips.length,
      clock: results.clockModel !== undefined && results.augurClock !== undefined,
      trace: (results.trace?.length ?? 0) > 0,
      skyline: results.skyline?.length,
      problems: results.problems,
    }).toStrictEqual({ tree: 86, clock: true, trace: true, skyline: 20, problems: [] });
  });

  test("each command loads its outputs without problems", async () => {
    const runs: ReadonlyArray<[AppCommand, string]> = [
      ["clock", "clock"],
      ["ancestral", "ancestral"],
      ["mugration", "mugration"],
      ["optimize", "optimize"],
      ["prune", "prune"],
    ];

    const loaded = await Promise.all(runs.map(async ([command, run]) => load(command, run)));

    expect(loaded.map((results) => [results.auspice?.tree.tips.length, results.problems])).toStrictEqual(
      runs.map(() => [86, []]),
    );
  });

  test("an unreadable output is reported with its path and the other outputs still load", async () => {
    const files: RunFileEntry[] = [
      { path: "timetree.auspice.json", size: 1, kind: "auspice" },
      { path: "timetree.clock-model.json", size: 1, kind: "clock-model" },
    ];

    const results = await loadRunResults("timetree", files, (path) =>
      Promise.resolve(path.endsWith(".auspice.json") ? fixtureText("timetree/timetree.auspice.json") : "{ not json"),
    );

    expect({
      tree: results.auspice !== undefined,
      problems: results.problems.map((problem) => problem.path),
    }).toStrictEqual({
      tree: true,
      problems: ["timetree.clock-model.json"],
    });
  });

  test("the zika mugration trait has more states than the palette and keeps the Auspice colours", () => {
    const text = fixtureText("mugration/mugration.auspice.json");

    expect(mutedAuspiceDocument(text, readAuspiceTree(parseAuspiceJson(text)))).toStrictEqual(parseAuspiceJson(text));
  });
});

function fixtureJson(path: string): unknown {
  return JSON.parse(fixtureText(path));
}

function fixtureTree(path: string): ResultTree {
  return readAuspiceTree(parseAuspiceJson(fixtureText(path)));
}

function datasetText(path: string): string {
  return readFileSync(new URL(path, DATA_ROOT), "utf8");
}

function byText(left: string, right: string): number {
  return left.localeCompare(right);
}

async function load(command: AppCommand, run: string) {
  const directory = new URL(`${run}/`, FIXTURE_ROOT);

  const files = readdirSync(directory).map((name): RunFileEntry => ({
    path: name,
    size: statSync(new URL(name, directory)).size,
    kind: KINDS.find(([pattern]) => pattern.test(name))?.[1] ?? null,
  }));

  return loadRunResults(command, files, (path) => Promise.resolve(fixtureText(`${run}/${path}`)));
}

function nodeMutations(nodeData: unknown): Map<string, string[]> {
  const { nodes } = z.object({ nodes: z.record(z.string(), z.object({ muts: z.array(z.string()) })) }).parse(nodeData);

  return new Map(Object.entries(nodes).map(([name, node]) => [name, node.muts.toSorted()]));
}

function rateVariance(clockModel: unknown): number {
  const model = z
    .object({ stats: z.object({ estimated: z.object({ cov: z.array(z.array(z.number())) }) }) })
    .parse(clockModel);

  const variance = model.stats.estimated.cov[0]?.[0];

  if (variance === undefined) {
    throw new Error("the clock model has no rate variance");
  }

  return variance;
}

function skylineJsonValues(json: unknown): number[][] {
  const coalescent = z
    .object({
      outputs: z.object({
        segments: z.array(
          z.object({ segment: z.object({ start: z.number(), end: z.number() }), N_e: z.object({ value: z.number() }) }),
        ),
      }),
    })
    .parse(json);

  return coalescent.outputs.segments.map((segment) => [segment.segment.start, segment.segment.end, segment.N_e.value]);
}

function fixtureText(path: string): string {
  return readFileSync(new URL(path, FIXTURE_ROOT), "utf8");
}
