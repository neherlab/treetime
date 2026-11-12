import { eventStream, FakeServer, json, runRecord } from "../../__tests__/api_server";
import { removeRun } from "../removeRun";

const CANCELLED = {
  seq: 1,
  time: "2026-09-25T10:00:02Z",
  type: "terminal",
  data: { status: "cancelled", job_id: "r1" },
};

describe("removal of one run", () => {
  test("a finished run is deleted with one request", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1": () => json(runRecord("r1", "ok")),
      "DELETE /api/runs/r1": () => new Response(null, { status: 204 }),
    });

    await removeRun(server.client(), "r1");

    expect(server.keys()).toStrictEqual(["GET /api/runs/r1", "DELETE /api/runs/r1"]);
  });

  test("a computing run is cancelled, and deleted once its terminal event arrives", async () => {
    const server = new FakeServer({
      "GET /api/runs/r1": () => json(runRecord("r1", "running")),
      "POST /api/runs/r1/cancel": () => json({ cancelled: true }),
      "GET /api/runs/r1/events?from=0": () => eventStream([CANCELLED]),
      "DELETE /api/runs/r1": () => new Response(null, { status: 204 }),
    });

    await removeRun(server.client(), "r1");

    expect(server.keys()).toStrictEqual([
      "GET /api/runs/r1",
      "POST /api/runs/r1/cancel",
      "GET /api/runs/r1/events?from=0",
      "DELETE /api/runs/r1",
    ]);
  });

  test("a run that no longer exists counts as removed", async () => {
    const server = new FakeServer({});

    await removeRun(server.client(), "r1");

    expect(server.keys()).toStrictEqual(["GET /api/runs/r1"]);
  });
});
