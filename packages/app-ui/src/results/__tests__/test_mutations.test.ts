import { describe, expect, test } from "vitest";

import {
  ancestorStates,
  branchesByMutationCount,
  recurrentSites,
  stateChanges,
  uncertainAncestors,
} from "../mutations";
import { parseAuspiceJson, readAuspiceTree } from "../tree";

const TREE = readAuspiceTree(
  parseAuspiceJson(
    JSON.stringify({
      meta: {},
      tree: {
        name: "root",
        node_attrs: { div: 0, country: { value: "a", confidence: { a: 0.6, b: 0.4 } } },
        children: [
          {
            name: "X",
            node_attrs: { div: 0, country: { value: "b", confidence: { b: 0.95 } } },
            branch_attrs: { mutations: { nuc: ["A10G", "C20T", "G30A"] } },
            children: [
              { name: "x1", node_attrs: { country: { value: "b" } }, branch_attrs: { mutations: { nuc: ["C20A"] } } },
              { name: "x2", node_attrs: { country: { value: "a" } } },
            ],
          },
          { name: "y", node_attrs: { country: { value: "b" } }, branch_attrs: { mutations: { nuc: ["A10T", "-5A"] } } },
        ],
      },
    }),
  ),
);

describe("mutation tables", () => {
  test("branches are ordered by their number of mutations", () => {
    expect(
      branchesByMutationCount(TREE).map((branch) => [branch.name, branch.tips, branch.mutations.length]),
    ).toStrictEqual([
      ["X", 2, 3],
      ["y", 1, 2],
      ["x1", 1, 1],
    ]);
  });

  test("sites mutated on more than one branch are listed with their branch count", () => {
    expect(recurrentSites(TREE)).toStrictEqual([
      { position: 10, branches: 2 },
      { position: 20, branches: 2 },
    ]);
  });

  test("deletions and insertions count at their position", () => {
    const withIndels = readAuspiceTree(
      parseAuspiceJson(
        JSON.stringify({
          meta: {},
          tree: {
            name: "root",
            children: [
              { name: "a", branch_attrs: { mutations: { nuc: ["A5-", "C7T"] } } },
              { name: "b", branch_attrs: { mutations: { nuc: ["-5G", "N7A"] } } },
            ],
          },
        }),
      ),
    );

    expect(recurrentSites(withIndels)).toStrictEqual([
      { position: 5, branches: 2 },
      { position: 7, branches: 2 },
    ]);
  });
});

describe("trait tables", () => {
  test("state changes count branches whose parent and child states differ", () => {
    expect(stateChanges(TREE, "country")).toStrictEqual([
      { from: "a", to: "b", branches: 2 },
      { from: "b", to: "a", branches: 1 },
    ]);
  });

  test("uncertain ancestors are those whose state has a probability below 0.8", () => {
    expect({
      all: ancestorStates(TREE, "country").map((state) => [state.name, state.probability]),
      uncertain: uncertainAncestors(TREE, "country").map((state) => state.name),
    }).toStrictEqual({
      all: [
        ["root", 0.6],
        ["X", 0.95],
      ],
      uncertain: ["root"],
    });
  });
});
