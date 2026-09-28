import type { StalePath } from "@neherlab/app-contracts";
import {
  cladeInRuns,
  runsEvents,
  runsFiles,
  runsGet,
  runsList,
  runsResults,
  StreamEndedError,
  version,
} from "@neherlab/app-contracts/client";
import { QueryClient, QueryObserver, type QueryKey } from "@tanstack/react-query";
import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import { followAppEvents, invalidateStale, runEventsQueryOptions, runEventStream } from "../api/events";
import { requestKey } from "../api/keys";
import { EMPTY_PROGRESS, foldRunEvents, type RunEvent } from "../results/progress";
import { eventStream, FakeServer, json, noDelay, RECORD } from "./api_server";

const TIME = "2026-09-25T10:00:01Z";

const RUN_EVENTS: RunEvent[] = [
  { seq: 0, time: TIME, type: "started", data: { job_id: "r1", command: "clock" } },
  { seq: 1, time: TIME, type: "progress", data: { stage: "clock", fraction: 0.5, message: "fitting" } },
  { seq: 2, time: TIME, type: "log", data: { level: "warn", message: "few dates" } },
  {
    seq: 3,
    time: TIME,
    type: "terminal",
    data: {
      status: "ok",
      job_id: "r1",
      result: { command: "clock", output_files: [{ path: "out/clock.nwk", kind: "nwk" }] },
    },
  },
];

const RUN_ABC_STALE: StalePath[] = [
  { path: "/api/runs", scope: "exact" },
  { path: "/api/runs/abc", scope: "subtree" },
  { path: "/api/clade-in-runs", scope: "exact" },
];

const RESYNC_STALE: StalePath[] = [
  { path: "/api/runs", scope: "subtree" },
  { path: "/api/clade-in-runs", scope: "exact" },
];

function appEvent(seq: number, kind: "resync" | "run-updated", stale: StalePath[]) {
  return kind === "resync"
    ? { seq, time: "2026-09-25T10:00:00Z", stale, kind }
    : {
        seq,
        time: "2026-09-25T10:00:00Z",
        stale,
        kind,
        run: {
          id: "abc",
          title: "abc",
          command: "clock",
          status: "ok",
          pinned: false,
          created_at: "2026-09-25T08:00:00Z",
          changed_settings: [],
          headline: {},
        },
      };
}

function seeded(keys: readonly QueryKey[]): QueryClient {
  const queryClient = new QueryClient();

  for (const key of keys) {
    queryClient.setQueryData(key, { cached: key });
  }

  return queryClient;
}

function invalidated(queryClient: QueryClient): string[] {
  return queryClient
    .getQueryCache()
    .getAll()
    .filter((query) => query.state.isInvalidated)
    .map((query) => JSON.stringify(query.queryKey))
    .toSorted();
}

describe("api events invalidation", () => {
  const client = new FakeServer({}).client();

  const runList = requestKey(client, runsList);
  const runListWithQuery = requestKey(client, ({ client }) => client.get({ url: "/api/runs", query: { page: 2 } }));
  const abc = requestKey(client, (context) => runsGet({ ...context, path: { id: "abc" } }));
  const abcResults = requestKey(client, (context) => runsResults({ ...context, path: { id: "abc" } }));
  const abcFiles = requestKey(client, (context) => runsFiles({ ...context, path: { id: "abc" } }));
  const abcd = requestKey(client, (context) => runsGet({ ...context, path: { id: "abcd" } }));
  const xyzResults = requestKey(client, (context) => runsResults({ ...context, path: { id: "xyz" } }));
  const clade = requestKey(client, (context) => cladeInRuns({ ...context, body: { run: "xyz", node: "NODE_1" } }));
  const versionKey = requestKey(client, version);

  const keys = [runList, runListWithQuery, abc, abcResults, abcFiles, abcd, xyzResults, clade, versionKey];

  test("a change of one run invalidates the run list, that run and its outputs, and the clade computation", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, RUN_ABC_STALE);

    expect(invalidated(queryClient)).toStrictEqual(
      [runList, runListWithQuery, abc, abcResults, abcFiles, clade].map((key) => JSON.stringify(key)).toSorted(),
    );
  });

  test("a change of one run leaves the outputs of other runs valid", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, RUN_ABC_STALE);

    expect([xyzResults, abcd, versionKey].map((key) => queryClient.getQueryState(key)?.isInvalidated)).toStrictEqual([
      false,
      false,
      false,
    ]);
  });

  test("a stale run subtree invalidates the run and its outputs and nothing else", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, [{ path: "/api/runs/abc", scope: "subtree" }]);

    expect(invalidated(queryClient)).toStrictEqual(
      [abc, abcResults, abcFiles].map((key) => JSON.stringify(key)).toSorted(),
    );
  });

  test("a resync invalidates every run query and the clade computation", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, RESYNC_STALE);

    expect(invalidated(queryClient)).toStrictEqual(
      keys.flatMap((key) => (key === versionKey ? [] : [JSON.stringify(key)])).toSorted(),
    );
  });

  test("an event without stale paths invalidates nothing", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, []);

    expect(invalidated(queryClient)).toStrictEqual([]);
  });
});

describe("api events app stream", () => {
  test("the app stream starts from the beginning, invalidates the stale paths of each event, and resumes after a drop", async () => {
    const controller = new AbortController();
    const xyz: StalePath[] = [{ path: "/api/runs/xyz", scope: "subtree" }];

    const server = new FakeServer({
      "GET /api/events?from=0": () =>
        eventStream(
          [appEvent(7, "resync", RESYNC_STALE), appEvent(8, "run-updated", RUN_ABC_STALE)],
          new TypeError("connection reset"),
        ),
      "GET /api/events?from=9": () => eventStream([appEvent(9, "run-updated", xyz)]),
    });

    const client = server.client();
    const runList = requestKey(client, runsList);
    const abcResults = requestKey(client, (context) => runsResults({ ...context, path: { id: "abc" } }));
    const xyzResults = requestKey(client, (context) => runsResults({ ...context, path: { id: "xyz" } }));
    const keys = [runList, abcResults, xyzResults];
    const queryClient = seeded(keys);
    const invalidate = queryClient.invalidateQueries.bind(queryClient);
    const rounds: string[][] = [];

    queryClient.invalidateQueries = (filters) => {
      const done = invalidate(filters);
      rounds.push(invalidated(queryClient));

      for (const key of keys) {
        queryClient.setQueryData(key, { cached: key });
      }

      if (rounds.length === 3) {
        controller.abort();
      }

      return done;
    };

    await followAppEvents({ client, queryClient, signal: controller.signal, sleep: noDelay });

    expect(server.keys()).toStrictEqual(["GET /api/events?from=0", "GET /api/events?from=9"]);
    expect(rounds).toStrictEqual([
      [runList, abcResults, xyzResults].map((key) => JSON.stringify(key)).toSorted(),
      [runList, abcResults].map((key) => JSON.stringify(key)).toSorted(),
      [JSON.stringify(xyzResults)],
    ]);
  });

  test("the app stream keeps reconnecting while the server is away", async () => {
    const controller = new AbortController();

    const server = new FakeServer({
      "GET /api/events?from=0": (_request, attempt) => {
        if (attempt === 20) {
          controller.abort();
        }

        return json({ code: "internal_error", message: "restarting", causes: [] }, 503);
      },
    });

    await followAppEvents({
      client: server.client(),
      queryClient: new QueryClient(),
      signal: controller.signal,
      sleep: noDelay,
    });

    expect(server.sent).toHaveLength(21);
  });

  test("an invalid app event fails the app stream", async () => {
    const server = new FakeServer({ "GET /api/events?from=0": () => eventStream([{ seq: 1 }]) });

    await expect(
      followAppEvents({
        client: server.client(),
        queryClient: new QueryClient(),
        signal: new AbortController().signal,
        sleep: noDelay,
      }),
    ).rejects.toBeInstanceOf(ZodError);
    expect(server.sent).toHaveLength(1);
  });
});

describe("api events run stream", () => {
  test("the run events query folds the stream into the run progress up to the terminal event", async () => {
    const server = new FakeServer({ "GET /api/runs/r1/events?from=0": () => eventStream(RUN_EVENTS) });
    const queryClient = new QueryClient();

    const progress = await queryClient.query(runEventsQueryOptions(server.client(), "r1", { sleep: noDelay }));

    expect(progress).toStrictEqual(foldRunEvents(EMPTY_PROGRESS, RUN_EVENTS));
    expect(progress.terminal).toStrictEqual(RUN_EVENTS[3]?.data);
    expect(queryClient.getQueryData(["api", "runs", "r1", "events"])).toStrictEqual(progress);
  });

  test("a dropped run stream resumes after the last event it delivered", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events?from=0": () => eventStream(RUN_EVENTS.slice(0, 2), new TypeError("connection reset")),
      "GET /api/runs/r1/events?from=2": () => eventStream(RUN_EVENTS.slice(2)),
    });

    const progress = await new QueryClient().query(runEventsQueryOptions(server.client(), "r1", { sleep: noDelay }));

    expect(progress).toStrictEqual(foldRunEvents(EMPTY_PROGRESS, RUN_EVENTS));
    expect(server.keys()).toStrictEqual(["GET /api/runs/r1/events?from=0", "GET /api/runs/r1/events?from=2"]);
  });

  test("a run stream that stops delivering events without a terminal event fails the query", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events?from=0": () => eventStream(RUN_EVENTS.slice(0, 2)),
      "GET /api/runs/r1/events?from=2": () => eventStream([]),
    });

    const error = await new QueryClient()
      .query(runEventsQueryOptions(server.client(), "r1", { sleep: noDelay }))
      .catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(StreamEndedError);
    expect(server.sent).toHaveLength(9);
  });

  test("aborting a run stream ends it quietly after the events it delivered", async () => {
    const controller = new AbortController();
    const server = new FakeServer({ "GET /api/runs/r1/events?from=0": () => eventStream(RUN_EVENTS) });
    const received: RunEvent[] = [];

    for await (const event of runEventStream(server.client(), "r1", { signal: controller.signal, sleep: noDelay })) {
      received.push(event);
      controller.abort();
    }

    expect(received).toStrictEqual(RUN_EVENTS.slice(0, 1));
  });

  test("a run stream of an unknown run fails with the ApiError without reconnecting", async () => {
    const server = new FakeServer({});

    const received: RunEvent[] = [];

    await expect(
      (async () => {
        for await (const event of runEventStream(server.client(), "r9", { sleep: noDelay })) {
          received.push(event);
        }
      })(),
    ).rejects.toMatchObject({ status: 404 });
    expect(received).toStrictEqual([]);
    expect(server.keys()).toStrictEqual(["GET /api/runs/r9/events?from=0"]);
  });

  test("invalidating the run does not restart its finished event stream", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1/events?from=0": () => eventStream(RUN_EVENTS),
      "GET /api/runs/r1": () => json(RECORD),
    });

    const client = server.client();
    const queryClient = new QueryClient();
    const options = runEventsQueryOptions(client, "r1", { sleep: noDelay });
    const observer = new QueryObserver(queryClient, options);
    const unsubscribe = observer.subscribe(() => undefined);

    await expect.poll(() => observer.getCurrentResult().data?.terminal).toBeDefined();
    await invalidateStale(queryClient, [{ path: "/api/runs/r1", scope: "subtree" }]);
    unsubscribe();

    expect(server.keys()).toStrictEqual(["GET /api/runs/r1/events?from=0"]);
    expect(requestKey(client, (context) => runsEvents({ ...context, path: { id: "r1" } }))).toStrictEqual(
      options.queryKey,
    );
  });
});
