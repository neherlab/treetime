import type { RunRecord, UiDraft } from "@neherlab/app-contracts";

import { useDraftStore } from "../draft";
import { heldUploadRuns, isPristine, ownDraft, releasedUploadRuns, runOptionsDraft } from "../draftRules";
import { freshDraft } from "../draftSchema";
import { blankDraftState } from "./draftState";

const FRESH = freshDraft("timetree");

const OWN: UiDraft = {
  ...FRESH,
  config: { ...FRESH.config, tree: "/runs/.inputs/a/tree.nwk", clock_rate: 0.001 },
  sources: { tree: { label: "tree.nwk", origin: { kind: "upload" }, size: 10, run_id: "u1" } },
  upload_run_id: "u1",
  search: "rate",
};

const RUN: RunRecord = {
  id: "r1",
  title: "Baseline",
  renamed: false,
  command: "timetree",
  config: { tree: "inputs/tree.nwk", clock_filter: 3, output_all: "out" },
  status: "ok",
  prepared: true,
  pinned: false,
  group_id: "r1",
  size: 30,
  seed: 4_242,
  created_at: "2026-11-12T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [{ setting: "tree", path: "inputs/tree.nwk", size: 10, sha256: "a" }],
  input_sets: [],
  uploads: [],
  changed_settings: ["clock_filter"],
  headline: {},
  output_files: [],
  warnings: [],
};

const RUN_OPTIONS = runOptionsDraft(RUN, { tree: "/runs/.inputs/a/tree.nwk", clock_filter: 3 }, OWN) ?? FRESH;

beforeEach(() => {
  useDraftStore.setState({ ...blankDraftState("timetree"), draft: OWN, epoch: 0 });
});

describe("options of a shown run", () => {
  test("hold the inputs and changed settings of the run, its seed, and the display of the own draft", () => {
    expect(RUN_OPTIONS).toStrictEqual({
      ...FRESH,
      config: { ...FRESH.config, tree: "/runs/.inputs/a/tree.nwk", clock_filter: 3, seed: 4_242 },
      sources: { tree: { label: "tree.nwk", origin: { kind: "run" }, run_id: "r1" } },
      from_run_id: "r1",
      search: "rate",
    });
  });

  test("come only from a run whose settings the run layer has prepared", () => {
    expect(runOptionsDraft({ ...RUN, prepared: false }, { clock_filter: 3 }, OWN)).toBeUndefined();
  });

  test("keep the example dataset of the run", () => {
    expect(runOptionsDraft({ ...RUN, example: "ebola/100" }, { clock_filter: 3 }, OWN)?.example).toBe("ebola/100");
  });

  test("showing a run parks the own draft, and leaving the run restores it unchanged", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    const viewed = useDraftStore.getState();
    useDraftStore.getState().leaveRun();
    const left = useDraftStore.getState();

    expect([viewed.draft, ownDraft(viewed), left.draft, left.viewing, left.undo]).toStrictEqual([
      RUN_OPTIONS,
      OWN,
      OWN,
      undefined,
      undefined,
    ]);
  });

  test("an edit that leaves the settings as they are keeps the run read-only", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config });

    expect(useDraftStore.getState().viewing?.runId).toBe("r1");
  });

  test("the first edit turns the options into a variant draft of the run with its seed and keeps the own draft", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    const { epoch } = useDraftStore.getState();
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config, coalescent: "skyline" });
    const state = useDraftStore.getState();

    expect([state.draft, state.viewing, state.epoch, state.undo, ownDraft(state)]).toStrictEqual([
      { ...RUN_OPTIONS, config: { ...RUN_OPTIONS.config, coalescent: "skyline" } },
      undefined,
      epoch,
      undefined,
      OWN,
    ]);
  });

  test("a reset of the variant shows the own draft again, and undo restores the variant", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().pickChoice({ choice: "root", option: "keep" });
    const variant = useDraftStore.getState().draft;
    useDraftStore.getState().discard("reset", "timetree");
    const reset = useDraftStore.getState();
    useDraftStore.getState().undoDiscard(reset.undo?.token ?? Symbol("missing"), undefined);
    const undone = useDraftStore.getState();

    expect([reset.draft, reset.parked, undone.draft, ownDraft(undone)]).toStrictEqual([OWN, undefined, variant, OWN]);
  });

  test("a reset of an unedited shown run shows the own draft again without an undo", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().discard("reset", "timetree");
    const state = useDraftStore.getState();

    expect([state.draft, state.viewing, state.undo]).toStrictEqual([OWN, undefined, undefined]);
  });

  test("a new analysis from a shown run with an untouched own draft offers no undo", () => {
    useDraftStore.setState({ draft: FRESH });
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().discard("new", "timetree");

    expect([useDraftStore.getState().draft.sources, useDraftStore.getState().undo]).toStrictEqual([{}, undefined]);
  });

  test("starting the variant shows the own draft again", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config, coalescent: "skyline" });
    useDraftStore.getState().finishStart();
    const state = useDraftStore.getState();

    expect([state.draft, state.parked]).toStrictEqual([OWN, undefined]);
  });

  test("opening another run while a variant is edited offers to undo the variant and keeps the own draft", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config, coalescent: "skyline" });
    useDraftStore.getState().viewRun("r2", { ...RUN_OPTIONS, from_run_id: "r2" });
    const state = useDraftStore.getState();

    expect([state.undo?.reason, ownDraft(state)]).toStrictEqual(["replaced", OWN]);
  });

  test("a variant of a pristine draft offers no undo", () => {
    useDraftStore.setState({ draft: FRESH });
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().setConfig({ ...RUN_OPTIONS.config, coalescent: "skyline" });

    expect(useDraftStore.getState().undo).toBeUndefined();
  });

  test("a draft that is already a variant of the shown run stays editable", () => {
    const variant = { ...RUN_OPTIONS, config: { ...RUN_OPTIONS.config, coalescent: "skyline" } };
    useDraftStore.setState({ draft: variant });
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);

    expect([useDraftStore.getState().draft, useDraftStore.getState().viewing]).toStrictEqual([variant, undefined]);
  });
});

describe("revert to the settings of a shown run", () => {
  test("shows the run read-only again, and undo restores the edited variant", () => {
    const variant = { ...RUN_OPTIONS, config: { ...RUN_OPTIONS.config, coalescent: "skyline" } };
    useDraftStore.setState({ draft: variant });
    useDraftStore.getState().revertToRun("r1", RUN_OPTIONS);
    const reverted = useDraftStore.getState();
    useDraftStore.getState().undoDiscard(reverted.undo?.token ?? Symbol("missing"), "r1");
    const undone = useDraftStore.getState();

    expect([
      reverted.draft,
      reverted.viewing?.runId,
      reverted.undo?.reason,
      undone.draft,
      undone.viewing,
    ]).toStrictEqual([RUN_OPTIONS, "r1", "revert", variant, undefined]);
  });
});

describe("reset with undo", () => {
  test("reset starts from the defaults with the display kept, and undo restores the draft", () => {
    useDraftStore.getState().discard("reset", "clock");
    const reset = useDraftStore.getState();
    useDraftStore.getState().undoDiscard(reset.undo?.token ?? Symbol("missing"), undefined);
    const undone = useDraftStore.getState();

    expect([reset.draft, reset.undo?.reason, undone.draft, undone.undo]).toStrictEqual([
      { ...freshDraft("clock"), search: "rate" },
      "reset",
      OWN,
      undefined,
    ]);
  });

  test("a reset of a pristine draft offers no undo", () => {
    useDraftStore.setState({ draft: FRESH });
    useDraftStore.getState().discard("reset", "timetree");

    expect(useDraftStore.getState().undo).toBeUndefined();
  });

  test("undo of a new analysis started from a shown run restores the own draft, without the run it no longer shows", () => {
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    useDraftStore.getState().discard("new", "timetree");
    useDraftStore.getState().undoDiscard(useDraftStore.getState().undo?.token ?? Symbol("missing"), undefined);
    const state = useDraftStore.getState();

    expect([state.draft, state.viewing]).toStrictEqual([OWN, undefined]);
  });

  test("an undo token of an earlier discard changes nothing", () => {
    useDraftStore.getState().discard("reset", "timetree");
    const stale = useDraftStore.getState().undo?.token ?? Symbol("missing");
    useDraftStore.getState().setConfig({ ...FRESH.config, clock_rate: 0.002 });
    useDraftStore.getState().discard("reset", "timetree");
    const before = useDraftStore.getState();
    useDraftStore.getState().undoDiscard(stale, undefined);
    const after = useDraftStore.getState();

    expect([after.draft, after.undo, stale === before.undo?.token]).toStrictEqual([before.draft, before.undo, false]);
  });
});

describe("uploads in flight while a run is shown", () => {
  test("opening a run keeps the file choices of the own draft, and a finished upload lands in it", () => {
    const selection = Symbol("pick");
    useDraftStore.getState().claimSelections(["alignment"], selection);
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    const alignment = { label: "aln.fasta", origin: { kind: "upload" as const }, size: 20, run_id: "u1" };

    useDraftStore.getState().assignInputs(
      {
        config: { ...OWN.config, alignment: ["/runs/.inputs/b/aln.fasta"] },
        sources: { ...OWN.sources, alignment },
      },
      ["alignment"],
    );

    const state = useDraftStore.getState();

    expect([state.draft, state.viewing?.runId, state.parked?.draft.sources["alignment"]]).toStrictEqual([
      RUN_OPTIONS,
      "r1",
      alignment,
    ]);
  });
});

describe("upload runs of the form", () => {
  test("a pristine draft holds no upload run, the own draft holds its run", () => {
    expect([isPristine(FRESH), isPristine(OWN), [...heldUploadRuns(useDraftStore.getState())]]).toStrictEqual([
      true,
      false,
      ["u1"],
    ]);
  });

  test("showing a run and leaving it release no upload run", () => {
    const before = useDraftStore.getState();
    useDraftStore.getState().viewRun("r1", RUN_OPTIONS);
    const viewed = useDraftStore.getState();
    useDraftStore.getState().leaveRun();

    expect([releasedUploadRuns(before, viewed), releasedUploadRuns(viewed, useDraftStore.getState())]).toStrictEqual([
      [],
      [],
    ]);
  });

  test("a discard releases the upload run only when its undo closes", () => {
    const before = useDraftStore.getState();
    useDraftStore.getState().discard("reset", "timetree");
    const discarded = useDraftStore.getState();
    useDraftStore.getState().dropUndo(discarded.undo?.token ?? Symbol("missing"));

    expect([
      releasedUploadRuns(before, discarded),
      releasedUploadRuns(discarded, useDraftStore.getState()),
    ]).toStrictEqual([[], ["u1"]]);
  });

  test("removing the last uploaded file releases the upload run", () => {
    const before = useDraftStore.getState();
    useDraftStore.getState().setSource("tree", undefined);
    const after = useDraftStore.getState();

    expect([after.draft.upload_run_id, releasedUploadRuns(before, after)]).toStrictEqual([undefined, ["u1"]]);
  });

  test("removing a file while another file is chosen keeps the upload run", () => {
    useDraftStore.getState().claimSelections(["alignment"], Symbol("pick"));
    useDraftStore.getState().setSource("tree", undefined);

    expect(useDraftStore.getState().draft.upload_run_id).toBe("u1");
  });

  test("replacing every upload with files of the computer releases the upload run", () => {
    const before = useDraftStore.getState();

    useDraftStore.getState().assignInputs(
      {
        config: { ...OWN.config, tree: "/data/tree.nwk" },
        sources: { tree: { label: "tree.nwk", origin: { kind: "local" } } },
      },
      ["tree"],
    );

    expect(releasedUploadRuns(before, useDraftStore.getState())).toStrictEqual(["u1"]);
  });
});
