import type { RunRecord } from "@neherlab/app-contracts";
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
import { QueryClient, partialMatchKey } from "@tanstack/react-query";
import { describe, expect, expectTypeOf, test } from "vitest";

import { apiKey, apiQueryOptions } from "../api/hooks";
import { pathKey, requestKey, type ApiRequest } from "../api/keys";
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

  test("a stale path is the key prefix of every request below it and of nothing else", () => {
    const prefix = pathKey("/api/runs/abc");

    expect(prefix).toStrictEqual(key((context) => runsGet({ ...context, path: { id: "abc" } })));
    expect(
      partialMatchKey(
        key((context) => runsResults({ ...context, path: { id: "abc" } })),
        prefix,
      ),
    ).toBe(true);
    expect(
      partialMatchKey(
        key((context) => runsCompare({ ...context, path: { id: "abc", other: "x" } })),
        prefix,
      ),
    ).toBe(true);
    expect(
      partialMatchKey(
        key((context) => runsGet({ ...context, path: { id: "abcd" } })),
        prefix,
      ),
    ).toBe(false);
    expect(partialMatchKey(key(runsList), prefix)).toBe(false);
  });

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
});
