import type { JobId } from "@neherlab/app-contracts";

import { hiddenRunIds, useRunDeletionStore } from "../runDeletion";

beforeEach(() => {
  useRunDeletionStore.setState({ deletions: [], removed: [] });
});

describe("deletion of runs with Undo", () => {
  test("hiding runs takes them out of the list at once and sends nothing", () => {
    const sent = recorder();

    useRunDeletionStore.getState().hide(["a", "b"]);

    expect([...hidden()]).toStrictEqual(["a", "b"]);
    expect(sent.ids).toStrictEqual([]);
  });

  test("undo shows exactly the hidden set again and a later commit sends no request", async () => {
    const sent = recorder();
    const { hide, undo, commit } = useRunDeletionStore.getState();
    const first = token(hide(["a", "b"]));
    hide(["c"]);

    undo(first);
    const outcome = await commit(first, sent.remove);

    expect([...hidden()]).toStrictEqual(["c"]);
    expect(outcome).toStrictEqual({ removed: [], failed: [] });
    expect(sent.ids).toStrictEqual([]);
  });

  test("closing the toast sends one request per run and keeps the runs hidden", async () => {
    const sent = recorder();
    const { hide, commit } = useRunDeletionStore.getState();
    const deletion = token(hide(["a", "b", "a"]));

    const outcome = await commit(deletion, sent.remove);
    await commit(deletion, sent.remove);

    expect(sent.ids).toStrictEqual(["a", "b"]);
    expect(outcome).toStrictEqual({ removed: ["a", "b"], failed: [] });
    expect([...hidden()]).toStrictEqual(["a", "b"]);
  });

  test("a run that is already hidden is not hidden or sent a second time", async () => {
    const sent = recorder();
    const { hide, commit } = useRunDeletionStore.getState();
    const first = token(hide(["a"]));
    const second = token(hide(["a", "b"]));

    await commit(first, sent.remove);
    await commit(second, sent.remove);

    expect(hide(["a", "b"])).toBeUndefined();
    expect(sent.ids).toStrictEqual(["a", "b"]);
  });

  test("a run whose delete fails shows again and is reported", async () => {
    const failure = new Error("run `b` is computing; cancel it before deleting it");
    const sent = recorder(new Map([["b", failure]]));
    const { hide, commit } = useRunDeletionStore.getState();

    const outcome = await commit(token(hide(["a", "b"])), sent.remove);

    expect(outcome).toStrictEqual({ removed: ["a"], failed: [{ id: "b", error: failure }] });
    expect([...hidden()]).toStrictEqual(["a"]);
  });

  test("the run list forgets a deleted run once the list no longer holds it", async () => {
    const sent = recorder();
    const { hide, commit, forget } = useRunDeletionStore.getState();
    await commit(token(hide(["a", "b"])), sent.remove);

    forget(new Set(["b", "c"]));

    expect([...hidden()]).toStrictEqual(["b"]);
  });
});

function hidden(): ReadonlySet<JobId> {
  return hiddenRunIds(useRunDeletionStore.getState());
}

function recorder(failures: ReadonlyMap<JobId, Error> = new Map()) {
  const ids: JobId[] = [];

  return {
    ids,
    remove: (id: JobId): Promise<void> => {
      ids.push(id);
      const failure = failures.get(id);

      return failure === undefined ? Promise.resolve() : Promise.reject(failure);
    },
  };
}

function token(hidden: symbol | undefined): symbol {
  if (hidden === undefined) {
    throw new Error("the runs were hidden already");
  }

  return hidden;
}
