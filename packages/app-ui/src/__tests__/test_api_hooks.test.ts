import { ApiError, runsCancel, runsGet, runsStart } from "@neherlab/app-contracts/client";
import { MutationObserver, QueryClient } from "@tanstack/react-query";
import { describe, expect, test } from "vitest";

import { apiKey, apiMutationOptions } from "../api/hooks";
import { FakeServer, json, RECORD } from "./api_server";

const NOT_FOUND = { code: "not_found", message: "When reading run `r9`", causes: ["no run with id `r9`"] };

describe("api mutations", () => {
  test("a mutation sends the request and seeds the cache entry its seed names with the response", async () => {
    const started = { ...RECORD, status: "running" };
    const server = new FakeServer({ "POST /api/runs/r1/start": () => json(started) });
    const client = server.client();
    const queryClient = new QueryClient();

    const mutation = new MutationObserver(
      queryClient,
      apiMutationOptions(
        client,
        queryClient,
        (context, { id, config }: { id: string; config: unknown }) =>
          runsStart({ ...context, path: { id }, body: { config } }),
        {
          seed:
            (_record, { id }) =>
            (context) =>
              runsGet({ ...context, path: { id } }),
        },
      ),
    );

    await expect(mutation.mutate({ id: "r1", config: { tree: "u.nwk" } })).resolves.toStrictEqual(started);
    expect(
      queryClient.getQueryData(apiKey(client, (context) => runsGet({ ...context, path: { id: "r1" } }))),
    ).toStrictEqual(started);
    expect(server.sent).toStrictEqual([{ key: "POST /api/runs/r1/start", body: '{"config":{"tree":"u.nwk"}}' }]);
  });

  test("a mutation without a seed leaves the cache unchanged", async () => {
    const server = new FakeServer({ "POST /api/runs/r1/cancel": () => json({ cancelled: true }) });
    const client = server.client();
    const queryClient = new QueryClient();

    const mutation = new MutationObserver(
      queryClient,
      apiMutationOptions(client, queryClient, (context, id: string) => runsCancel({ ...context, path: { id } })),
    );

    await expect(mutation.mutate("r1")).resolves.toStrictEqual({ cancelled: true });

    expect(queryClient.getQueryCache().getAll()).toStrictEqual([]);
  });

  test("a failed mutation rejects with the ApiError and seeds nothing", async () => {
    const server = new FakeServer({ "POST /api/runs/r9/start": () => json(NOT_FOUND, 404) });
    const client = server.client();
    const queryClient = new QueryClient();

    const mutation = new MutationObserver(
      queryClient,
      apiMutationOptions(
        client,
        queryClient,
        (context, id: string) => runsStart({ ...context, path: { id }, body: {} }),
        {
          seed: (_record, id) => (context) => runsGet({ ...context, path: { id } }),
        },
      ),
    );

    const error = await mutation.mutate("r9").catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 404, response: NOT_FOUND });
    expect(queryClient.getQueryCache().getAll()).toStrictEqual([]);
  });
});
