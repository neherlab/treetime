import type { UiDraft } from "@neherlab/app-contracts";
import pWaitFor from "p-wait-for";

import { blankDraftState } from "../../store/__tests__/draftState";
import { useDraftStore } from "../../store/draft";
import { runOptionsDraft } from "../../store/draftRules";
import { freshDraft } from "../../store/draftSchema";
import { UploadRunJanitor } from "../uploadRunJanitor";

const FRESH = freshDraft("timetree");

const OWN: UiDraft = {
  ...FRESH,
  config: { ...FRESH.config, tree: "/runs/.inputs/a/tree.nwk" },
  sources: { tree: { label: "tree.nwk", origin: { kind: "upload" }, size: 10, run_id: "u1" } },
  upload_run_id: "u1",
};

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

  test("undo of a reset deletes nothing", () => {
    const retired = watch();
    useDraftStore.getState().discard("reset", "timetree");
    useDraftStore.getState().undoDiscard(useDraftStore.getState().undo?.token ?? Symbol("missing"), undefined);

    expect(retired.ids).toStrictEqual([]);
  });

  test("replacing every uploaded file deletes the upload run", async () => {
    const retired = watch();
    useDraftStore.getState().setSource("tree", { label: "tree.nwk", origin: { kind: "local" } });
    await pWaitFor(() => retired.ids.length > 0);

    expect(retired.ids).toStrictEqual(["u1"]);
  });

  test("leaving the form for a run and coming back deletes nothing", () => {
    const retired = watch();

    const options = runOptionsDraft(
      {
        id: "r1",
        title: "Baseline",
        command: "timetree",
        config: {},
        status: "ok",
        pinned: false,
        group_id: "r1",
        size: 0,
        created_at: "2026-11-12T08:00:00Z",
        treetime_version: "1.0.0",
        inputs: [],
        uploads: [],
        changed_settings: [],
        headline: {},
        output_files: [],
        warnings: [],
      },
      {},
      OWN,
    );

    useDraftStore.getState().viewRun("r1", options);

    useDraftStore.getState().leaveRun();

    expect(retired.ids).toStrictEqual([]);
  });
});

function watch() {
  const ids: string[] = [];

  janitors.push(
    new UploadRunJanitor(
      useDraftStore,
      (runId) => {
        ids.push(runId);

        return Promise.resolve();
      },
      (message) => {
        throw new Error(message);
      },
    ),
  );

  return { ids };
}
