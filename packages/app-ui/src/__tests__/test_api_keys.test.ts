import type { RunRecord, StalePath } from "@neherlab/app-contracts";
import {
  cladeInRuns,
  configCheck,
  events,
  runsCompare,
  runsEvents,
  runsFile,
  runsGet,
  runsList,
  runsResults,
  version,
  type ApiClient,
} from "@neherlab/app-contracts/client";
import { QueryClient } from "@tanstack/react-query";
import { describe, expect, expectTypeOf, test } from "vitest";

import { apiKey, apiQueryOptions, resetApiQueries } from "../api/hooks";
import { requestKey, requestUrl, staleCoversKey, type ApiRequest } from "../api/keys";
import { FakeServer, json, RECORD } from "./api_server";

const CLIENT = new FakeServer({}).client();

const BODY = { command: "clock" as const, text: "tree: t.nwk", input_facts: null };

function key(request: ApiRequest): readonly unknown[] {
  return requestKey(CLIENT, request);
}

describe("api keys", () => {
  test.each([
    { name: "a collection", request: (context) => runsList(context), expected: ["api", "runs"] },
    {
      name: "a resource",
      request: (context) => runsGet({ ...context, path: { id: "abc" } }),
      expected: ["api", "runs", "abc"],
    },
    {
      name: "a nested resource",
      request: (context) => runsResults({ ...context, path: { id: "abc" } }),
      expected: ["api", "runs", "abc", "results"],
    },
    {
      name: "two path parameters",
      request: (context) => runsCompare({ ...context, path: { id: "abc", other: "xyz" } }),
      expected: ["api", "runs", "abc", "compare", "xyz"],
    },
    {
      name: "an escaped path parameter",
      request: (context) => runsGet({ ...context, path: { id: "a b/c" } }),
      expected: ["api", "runs", "a b/c"],
    },
    { name: "an operation without parameters", request: (context) => version(context), expected: ["api", "version"] },
    {
      name: "a query parameter",
      request: (context) => runsFile({ ...context, path: { id: "abc" }, query: { path: "out/a b.nwk" } }),
      expected: ["api", "runs", "abc", "file", { path: "out/a b.nwk" }],
    },
    {
      name: "an event stream",
      request: (context) => runsEvents({ ...context, path: { id: "abc" } }),
      expected: ["api", "runs", "abc", "events"],
    },
    {
      name: "an event stream with a query",
      request: (context) => events({ ...context, query: { from: 3 } }),
      expected: ["api", "events", { from: 3 }],
    },
    {
      name: "a computation",
      request: (context) => configCheck({ ...context, body: BODY }),
      expected: ["api", "check-config", BODY],
    },
    {
      name: "a computation with a nested body",
      request: (context) => cladeInRuns({ ...context, body: { run: "abc", node: "NODE_1" } }),
      expected: ["api", "clade-in-runs", { run: "abc", node: "NODE_1" }],
    },
  ] satisfies Array<{ name: string; request: ApiRequest; expected: unknown[] }>)(
    "the key of $name is its path segments, then its query, then its body",
    ({ request, expected }) => {
      expect(key(request)).toStrictEqual(expected);
    },
  );

  test("query parameters are sorted by name and undefined parameters are dropped", () => {
    const listThings: ApiRequest = ({ client }) =>
      client.get({ url: "/api/things", query: { b: 2, c: undefined, a: "x" } });

    expect(key(listThings)).toStrictEqual(["api", "things", { a: "x", b: 2 }]);
    expect(JSON.stringify(key(listThings))).toBe('["api","things",{"a":"x","b":2}]');
  });

  test("deriving a key sends no request", () => {
    const server = new FakeServer({});

    requestKey(server.client(), (context) => runsGet({ ...context, path: { id: "abc" } }));
    requestKey(server.client(), (context) => runsEvents({ ...context, path: { id: "abc" } }));

    expect(server.sent).toStrictEqual([]);
  });

  test("a request that sends no request or two has no key", () => {
    expect(() => key(() => Promise.resolve({}))).toThrow("the request sent 0 requests to its client instead of one");
    expect(() =>
      key((context) => {
        void version(context);

        return version(context);
      }),
    ).toThrow("the request sent 2 requests to its client instead of one");
  });

  test.each([
    { name: "the run list", scope: "exact", path: "/api/runs", request: runsList, expected: true },
    {
      name: "the run list with a query",
      scope: "exact",
      path: "/api/runs",
      request: ({ client }) => client.get({ url: "/api/runs", query: { page: 2 } }),
      expected: true,
    },
    {
      name: "a run",
      scope: "exact",
      path: "/api/runs",
      request: (context) => runsGet({ ...context, path: { id: "abc" } }),
      expected: false,
    },
    {
      name: "a clade computation with a body",
      scope: "exact",
      path: "/api/clade-in-runs",
      request: (context) => cladeInRuns({ ...context, body: { run: "abc", node: "NODE_1" } }),
      expected: true,
    },
    {
      name: "a run's results",
      scope: "exact",
      path: "/api/runs/abc",
      request: (context) => runsResults({ ...context, path: { id: "abc" } }),
      expected: false,
    },
    {
      name: "a run file with a query",
      scope: "exact",
      path: "/api/runs/abc/file",
      request: (context) => runsFile({ ...context, path: { id: "abc" }, query: { path: "out/a.nwk" } }),
      expected: true,
    },
    {
      name: "the run",
      scope: "subtree",
      path: "/api/runs/abc",
      request: (context) => runsGet({ ...context, path: { id: "abc" } }),
      expected: true,
    },
    {
      name: "the run's results",
      scope: "subtree",
      path: "/api/runs/abc",
      request: (context) => runsResults({ ...context, path: { id: "abc" } }),
      expected: true,
    },
    {
      name: "a comparison of the run",
      scope: "subtree",
      path: "/api/runs/abc",
      request: (context) => runsCompare({ ...context, path: { id: "abc", other: "x" } }),
      expected: true,
    },
    {
      name: "a run whose id extends the run's id",
      scope: "subtree",
      path: "/api/runs/abc",
      request: (context) => runsGet({ ...context, path: { id: "abcd" } }),
      expected: false,
    },
    { name: "the run list", scope: "subtree", path: "/api/runs/abc", request: runsList, expected: false },
    {
      name: "another run's results",
      scope: "subtree",
      path: "/api/runs",
      request: (context) => runsResults({ ...context, path: { id: "xyz" } }),
      expected: true,
    },
    { name: "the version", scope: "subtree", path: "/api/runs", request: version, expected: false },
  ] satisfies Array<StalePath & { name: string; request: ApiRequest; expected: boolean }>)(
    "the $scope stale path $path covers $name: $expected",
    ({ scope, path, request, expected }) => {
      expect(staleCoversKey({ path, scope }, key(request))).toBe(expected);
    },
  );

  test("a query result is typed and cached under the key of its request", async () => {
    const server = new FakeServer({ "GET /api/runs/r1": () => json(RECORD) });
    const client: ApiClient = server.client();
    const queryClient = new QueryClient();
    const get = (context: Parameters<ApiRequest>[0]) => runsGet({ ...context, path: { id: "r1" } });

    const record = await queryClient.query(apiQueryOptions(client, get));

    expectTypeOf(record).toEqualTypeOf<RunRecord>();
    expect(record).toStrictEqual(RECORD);
    expect(queryClient.getQueryData(apiKey(client, get))).toStrictEqual(RECORD);
    expect(server.keys()).toStrictEqual(["GET /api/runs/r1"]);
  });

  test("the URL of a request carries its base, path and query", () => {
    const url = requestUrl(CLIENT, (context) =>
      runsFile({ ...context, path: { id: "r 1" }, query: { path: "out/a b.nwk" } }),
    );

    expect(url).toBe("http://treetime.test/api/runs/r%201/file?path=out%2Fa%20b.nwk");
  });
});

describe("api reset", () => {
  test("resetting the API queries keeps the queries of other sources", async () => {
    const queryClient = new QueryClient();
    const runs = key((context) => runsList(context));

    queryClient.setQueryData(runs, { runs: [], active_runs: 0 });
    queryClient.setQueryData(["preferences"], true);
    await resetApiQueries(queryClient);

    expect([queryClient.getQueryData(runs), queryClient.getQueryData(["preferences"])]).toStrictEqual([
      undefined,
      true,
    ]);
  });
});
