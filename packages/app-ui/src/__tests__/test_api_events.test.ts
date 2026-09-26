import {
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

function appEvent(seq: number, kind: "resync" | "run-deleted", stale: string[]) {
  return kind === "resync"
    ? { seq, time: "2026-09-25T10:00:00Z", stale, kind }
    : { seq, time: "2026-09-25T10:00:00Z", stale, kind, id: "abc" };
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

  const keys = [
    requestKey(client, runsList),
    requestKey(client, (context) => runsGet({ ...context, path: { id: "abc" } })),
    requestKey(client, (context) => runsResults({ ...context, path: { id: "abc" } })),
    requestKey(client, (context) => runsFiles({ ...context, path: { id: "abc" } })),
    requestKey(client, (context) => runsGet({ ...context, path: { id: "abcd" } })),
    requestKey(client, (context) => runsResults({ ...context, path: { id: "xyz" } })),
    requestKey(client, version),
  ];

  test("a stale run path invalidates the run and its outputs and nothing else", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, ["/api/runs/abc"]);

    expect(invalidated(queryClient)).toStrictEqual(
      ['["api","runs","abc"]', '["api","runs","abc","files"]', '["api","runs","abc","results"]'].toSorted(),
    );
  });

  test("the stale collection path invalidates everything below it", async () => {
    const queryClient = seeded(keys);

    await invalidateStale(queryClient, ["/api/runs"]);

    expect(invalidated(queryClient)).toStrictEqual(
      keys.flatMap((key) => (key[1] === "runs" ? [JSON.stringify(key)] : [])).toSorted(),
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

    const server = new FakeServer({
      "GET /api/events?from=0": () =>
        eventStream(
          [appEvent(7, "resync", ["/api/runs", "/api/clade-in-runs"]), appEvent(8, "run-deleted", ["/api/runs/abc"])],
          new TypeError("connection reset"),
        ),
      "GET /api/events?from=9": () => eventStream([appEvent(9, "run-deleted", ["/api/runs/xyz"])]),
    });

    const queryClient = new QueryClient();
    const stale: unknown[] = [];

    queryClient.invalidateQueries = (filters) => {
      stale.push(filters?.queryKey);

      if (stale.length === 4) {
        controller.abort();
      }

      return Promise.resolve();
    };

    await followAppEvents({ client: server.client(), queryClient, signal: controller.signal, sleep: noDelay });

    expect(server.keys()).toStrictEqual(["GET /api/events?from=0", "GET /api/events?from=9"]);
    expect(stale).toStrictEqual([
      ["api", "runs"],
      ["api", "clade-in-runs"],
      ["api", "runs", "abc"],
      ["api", "runs", "xyz"],
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
    await invalidateStale(queryClient, ["/api/runs/r1"]);
    unsubscribe();

    expect(server.keys()).toStrictEqual(["GET /api/runs/r1/events?from=0"]);
    expect(requestKey(client, (context) => runsEvents({ ...context, path: { id: "r1" } }))).toStrictEqual(
      options.queryKey,
    );
  });
});
