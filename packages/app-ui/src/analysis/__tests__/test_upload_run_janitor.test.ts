import type { UiDraft } from "@neherlab/app-contracts";
import pWaitFor from "p-wait-for";

import { blankDraftState } from "../../store/__tests__/draftState";
import { useDraftStore } from "../../store/draft";
import { runOptionsDraft } from "../../store/draftRules";
import { freshDraft } from "../../store/draftSchema";
import { type JanitorProblem, UploadRunJanitor } from "../uploadRunJanitor";

const FRESH = freshDraft("timetree");

const OWN: UiDraft = {
  ...FRESH,
  config: { ...FRESH.config, tree: "/runs/.inputs/a/tree.nwk" },
  sources: { tree: { label: "tree.nwk", origin: { kind: "upload" }, size: 10, run_id: "u1" } },
  upload_run_id: "u1",
};

const RUN_OPTIONS =
  runOptionsDraft(
    {
      id: "r1",
      title: "Baseline",
      renamed: false,
      command: "timetree",
      config: {},
      status: "ok",
      prepared: true,
      pinned: false,
      group_id: "r1",
      size: 0,
      created_at: "2026-11-12T08:00:00Z",
      treetime_version: "1.0.0",
      inputs: [],
      input_sets: [],
      uploads: [],
      changed_settings: [],
      headline: {},
      output_files: [],
      warnings: [],
    },
    {},
    OWN,
  ) ?? FRESH;

const janitors: UploadRunJanitor[] = [];

beforeEach(() => {
  useDraftStore.setState({ ...blankDraftState("timetree"), draft: OWN });
});

afterEach(() => {
  for (const janitor of janitors.splice(0)) {
    janitor.dispose();
  }
});

describe("upload run of a discarded form", () => {
  test("a reset deletes the upload run once, after its undo closes and its uploads end", async () => {
    const retired = watch();
    useDraftStore.getState().beginRunWork("u1");
    useDraftStore.getState().discard("reset", "timetree");
    const token = useDraftStore.getState().undo?.token ?? Symbol("missing");
    useDraftStore.getState().dropUndo(token);
    useDraftStore.getState().dropUndo(token);
    const whileUploading = [...retired.ids];
    useDraftStore.getState().endRunWork("u1");
    await pWaitFor(() => retired.ids.length > 0);

    expect([whileUploading, retired.ids]).toStrictEqual([[], ["u1"]]);
  });

  test("a reset while the upload run is being created deletes the run once its creation returns", async () => {
    useDraftStore.setState(blankDraftState("timetree"));
    const retired = watch();
    const { promise, resolve } = Promise.withResolvers<string>();

    const creation: Promise<string> = promise.then((id) => {
      useDraftStore.getState().settleUploadRun(creation, id);

      return id;
    });

    useDraftStore.getState().beginUploadRun(creation, undefined);
    useDraftStore.getState().discard("reset", "timetree");
    resolve("u2");
    await pWaitFor(() => retired.ids.length > 0);

    expect([retired.ids, useDraftStore.getState().draft.upload_run_id]).toStrictEqual([["u2"], undefined]);
  });

  test("undo of a reset deletes nothing", async () => {
    const retired = watch();
    useDraftStore.getState().discard("reset", "timetree");
    useDraftStore.getState().undoDiscard(useDraftStore.getState().undo?.token ?? Symbol("missing"), undefined);
    await retired.janitor.settled();

    expect(retired.ids).toStrictEqual([]);
  });

  test("a failed creation of the upload run of a discarded form is reported", async () => {
    useDraftStore.setState(blankDraftState("timetree"));
    const retired = watch();
    const { promise, reject } = Promise.withResolvers<string>();
    const creation = promise.finally(() => useDraftStore.getState().settleUploadRun(creation, undefined));
    useDraftStore.getState().beginUploadRun(creation, undefined);
    useDraftStore.getState().discard("reset", "timetree");
    reject(new Error("the back end stopped"));
    await retired.janitor.settled();

    expect([retired.ids, retired.problems]).toStrictEqual([
      [],
      [{ title: "The upload run of a discarded form was not created", description: "the back end stopped" }],
    ]);
  });

  test("opening a run while the upload run is being created keeps the run for the own draft", async () => {
    useDraftStore.setState(blankDraftState("timetree"));
    const retired = watch();
    const { promise, resolve } = Promise.withResolvers<string>();

    const creation: Promise<string> = promise.then((id) => {
      useDraftStore.getState().settleUploadRun(creation, id);

      return id;
    });

    useDraftStore.getState().beginUploadRun(creation, undefined);
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    resolve("u2");
    await creation;
    await retired.janitor.settled();
    useDraftStore.getState().leaveRun();

    expect([retired.ids, useDraftStore.getState().draft.upload_run_id]).toStrictEqual([[], "u2"]);
  });

  test("replacing every uploaded file deletes the upload run", async () => {
    const retired = watch();
    useDraftStore.getState().setSource("tree", { label: "tree.nwk", origin: { kind: "local" } });
    await pWaitFor(() => retired.ids.length > 0);

    expect(retired.ids).toStrictEqual(["u1"]);
  });

  test("leaving the form for a run and coming back deletes nothing", async () => {
    const retired = watch();
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().leaveRun();
    await retired.janitor.settled();

    expect(retired.ids).toStrictEqual([]);
  });

  test("an edit of a shown run keeps the upload run of the own draft", async () => {
    const retired = watch();
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config, clock_filter: 2 });
    await retired.janitor.settled();

    expect(retired.ids).toStrictEqual([]);
  });
});

function watch() {
  const ids: string[] = [];
  const problems: JanitorProblem[] = [];

  const janitor = new UploadRunJanitor(
    useDraftStore,
    (runId) => {
      ids.push(runId);

      return Promise.resolve();
    },
    (problem) => {
      problems.push(problem);
    },
  );

  janitors.push(janitor);

  return { ids, problems, janitor };
}
