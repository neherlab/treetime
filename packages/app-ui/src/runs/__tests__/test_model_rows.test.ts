import type { RunModels } from "@neherlab/app-contracts";

import { gammaText, indelText, modelColumnsInUse, modelRows, parameterSharing, type ModelRow } from "../modelRows";

const STATES = ["A", "C", "G", "T"];

const PI = [0.25, 0.25, 0.25, 0.25];

const RATES = [
  [0, 1, 1, 1],
  [1, 0, 1, 1],
  [1, 1, 0, 1],
  [1, 1, 1, 0],
];

const GENES: RunModels = {
  alignments: [{ alignment: 1, path: "/data/mpox/aln.fasta.xz", sites: 2000, indel_rate: 0.5 }],
  models: [
    {
      kind: "sequence",
      partition: 1,
      name: "OPG001",
      alignment: 1,
      path: "/data/mpox/aln.fasta.xz",
      sites: 1200,
      model: "infer",
      mu: 0.001_234_56,
      states: STATES,
      pi: PI,
      rates: RATES,
      // oxlint-disable-next-line anti-slop/no-shape-in-symbol-names -- `shape` is the gamma shape alpha of the generated contract
      gamma: { categories: 4, shape: 0.512_34, rates: [0.1, 0.5, 1, 2.4] },
      shared_params: ["exchange-rates", "gamma-shape"],
    },
    {
      kind: "sequence",
      partition: 2,
      alignment: 1,
      sites: 800,
      model: "hky85",
      mu: 0.002,
      states: STATES,
      pi: PI,
      rates: RATES,
      shared_params: [],
    },
    {
      kind: "sequence",
      partition: 3,
      name: "OPG003",
      alignment: 1,
      path: "/data/mpox/aln.fasta.xz",
      sites: 300,
      model: "infer",
      mu: 0.003,
      states: STATES,
      pi: PI,
      rates: RATES,
      shared_params: ["exchange-rates", "gamma-shape"],
    },
    {
      kind: "trait",
      attribute: "country",
      model: "infer",
      mu: 0.5,
      states: ["UK", "US"],
      pi: [0.5, 0.5],
      rates: RATES,
    },
  ],
  shared_params: ["exchange-rates", "gamma-shape"],
};

describe("model rows", () => {
  test("give each partition and trait its facts in document order", () => {
    const expected = [
      {
        id: "sequence-1",
        order: 0,
        label: "Partition 1 (OPG001)",
        alignment: "aln.fasta.xz",
        sites: 1200,
        model: "infer",
        mu: 0.001_234_56,
        indel: { perSite: 0.000_25, perAlignment: 0.5 },
        gamma: { alpha: 0.512_34, categories: 4 },
        shared: ["exchange-rates", "gamma-shape"],
      },
      {
        id: "sequence-2",
        order: 1,
        label: "Partition 2",
        alignment: "Alignment 1",
        sites: 800,
        model: "hky85",
        mu: 0.002,
        indel: { perSite: 0.000_25, perAlignment: 0.5 },
        gamma: undefined,
        shared: [],
      },
      {
        id: "sequence-3",
        order: 2,
        label: "Partition 3 (OPG003)",
        alignment: "aln.fasta.xz",
        sites: 300,
        model: "infer",
        mu: 0.003,
        indel: { perSite: 0.000_25, perAlignment: 0.5 },
        gamma: undefined,
        shared: ["exchange-rates", "gamma-shape"],
      },
      {
        id: "trait-country",
        order: 3,
        label: "Trait country",
        alignment: undefined,
        sites: undefined,
        model: "infer",
        mu: 0.5,
        indel: undefined,
        gamma: undefined,
        shared: [],
      },
    ];

    expect(modelRows(GENES).map(withoutSource)).toStrictEqual(expected);
  });

  test("leave the indel rate out when the run has no indel term", () => {
    const models: RunModels = {
      ...GENES,
      alignments: [{ alignment: 1, path: "/data/mpox/aln.fasta.xz", sites: 2000 }],
    };

    expect(modelRows(models).map((row) => row.indel)).toStrictEqual([undefined, undefined, undefined, undefined]);
  });

  test("give the indel rate per alignment only for an alignment without sites", () => {
    const models: RunModels = {
      ...GENES,
      alignments: [{ alignment: 1, path: "/data/mpox/aln.fasta.xz", sites: 0, indel_rate: 0.5 }],
    };

    expect(modelRows(models)[0]?.indel).toStrictEqual({ perSite: undefined, perAlignment: 0.5 });
  });
});

describe("model columns in use", () => {
  test("show the optional columns that some row fills", () => {
    expect(modelColumnsInUse(modelRows(GENES))).toStrictEqual({
      alignment: true,
      sites: true,
      indel: true,
      gamma: true,
      shared: true,
    });
  });

  test("hide the sequence columns of a run with traits only", () => {
    const traits: RunModels = { alignments: [], models: GENES.models.slice(3), shared_params: [] };

    expect(modelColumnsInUse(modelRows(traits))).toStrictEqual({
      alignment: false,
      sites: false,
      indel: false,
      gamma: false,
      shared: false,
    });
  });
});

describe("parameter sharing", () => {
  test("names the shared parameters and every partition that shares them", () => {
    const rows = modelRows(GENES);

    const expected = {
      params: "exchange rates, gamma shape",
      partitions: ["Partition 1 (OPG001)", "Partition 3 (OPG003)"],
      all: false,
    };

    expect(parameterSharing(rows, row(rows, "sequence-3"))).toStrictEqual(expected);
  });

  test("says that all partitions share when no partition has a named model", () => {
    const shared: RunModels = {
      ...GENES,
      models: GENES.models.filter((model) => model.kind === "trait" || model.shared_params.length > 0),
    };

    const rows = modelRows(shared);

    expect(parameterSharing(rows, row(rows, "sequence-1"))?.all).toBe(true);
  });

  test("is absent for a partition with a named model", () => {
    const rows = modelRows(GENES);

    expect(parameterSharing(rows, row(rows, "sequence-2"))).toBeUndefined();
  });
});

describe("model texts", () => {
  test.each([
    [{ perSite: 0.000_25, perAlignment: 0.5 }, "0.000250 per site (0.5000 per alignment)"],
    [{ perSite: undefined, perAlignment: 0.5 }, "0.5000 per alignment"],
  ])("writes the indel rate %o as %s", (indel, expected) => {
    expect(indelText(indel)).toBe(expected);
  });

  test("writes the gamma shape with the number of categories", () => {
    expect(gammaText({ alpha: 0.512_34, categories: 4 })).toBe("0.512 (K = 4)");
  });
});

function withoutSource({ source: _source, ...rest }: ModelRow) {
  return rest;
}

function row(rows: readonly ModelRow[], id: string): ModelRow {
  const found = rows.find((candidate) => candidate.id === id);

  if (found === undefined) {
    throw new Error(`no row ${id}`);
  }

  return found;
}
