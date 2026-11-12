import { freshDraft } from "../../store/draftSchema";
import { formUploads, leaveLosesWork, useTabWorkStore } from "../leaveGuard";

const NONE = { runs: 0, imports: 0, formUploads: 0 };

beforeEach(() => {
  useTabWorkStore.setState({ imports: 0, holdsWork: false });
});

describe("leave guard", () => {
  test.each([
    { name: "nothing to lose", work: NONE, expected: false },
    { name: "the tab holds runs", work: { ...NONE, runs: 2 }, expected: true },
    { name: "a package import is open", work: { ...NONE, imports: 1 }, expected: true },
    { name: "the form holds uploaded files", work: { ...NONE, formUploads: 1 }, expected: true },
  ])("$name: $expected", ({ work, expected }) => {
    expect(leaveLosesWork(work)).toBe(expected);
  });

  test("only uploaded files of the form count as files that no run uses", () => {
    const draft = freshDraft("clock");

    expect(
      formUploads({
        ...draft,
        sources: {
          tree: { label: "tree.nwk", origin: { kind: "upload" }, run_id: "r1" },
          metadata: { label: "dates.tsv", origin: { kind: "run" }, run_id: "r2" },
          alignment: { label: "aln.fasta", origin: { kind: "upload" }, run_id: "r1" },
        },
      }),
    ).toBe(2);
  });

  test("an import counts until it ends, and ending it twice counts once", () => {
    const { startImport } = useTabWorkStore.getState();
    const first = startImport();
    const second = startImport();

    first();
    first();
    const during = useTabWorkStore.getState().imports;
    second();

    expect([during, useTabWorkStore.getState().imports]).toStrictEqual([1, 0]);
  });
});
