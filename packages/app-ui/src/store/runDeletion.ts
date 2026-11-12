import type { JobId } from "@neherlab/app-contracts";
import { useMemo } from "react";
import { create } from "zustand";

interface Deletion {
  token: symbol;
  ids: readonly JobId[];
  sent: boolean;
}

interface RunDeletionState {
  deletions: readonly Deletion[];
  removed: readonly JobId[];
}

interface RunDeletionActions {
  hide: (ids: readonly JobId[]) => symbol | undefined;
  undo: (token: symbol) => void;
  commit: (token: symbol, remove: (id: JobId) => Promise<void>) => Promise<DeletionOutcome>;
  forget: (listed: ReadonlySet<JobId>) => void;
}

interface DeletionOutcome {
  removed: JobId[];
  failed: DeletionFailure[];
}

interface DeletionFailure {
  id: JobId;
  error: unknown;
}

type RunDeletionStore = RunDeletionState & RunDeletionActions;

export const useRunDeletionStore = create<RunDeletionStore>()((set, get) => ({
  deletions: [],
  removed: [],
  hide: (ids) => {
    const hidden = hiddenRunIds(get());
    const fresh = [...new Set(ids)].filter((id) => !hidden.has(id));

    if (fresh.length === 0) {
      return undefined;
    }

    const token = Symbol("delete");
    set(({ deletions }) => ({ deletions: [...deletions, { token, ids: fresh, sent: false }] }));

    return token;
  },
  undo: (token) => {
    set(({ deletions }) => ({
      deletions: deletions.filter((deletion) => deletion.token !== token || deletion.sent),
    }));
  },
  commit: async (token, remove) => {
    const deletion = get().deletions.find((candidate) => candidate.token === token && !candidate.sent);

    if (deletion === undefined) {
      return { removed: [], failed: [] };
    }

    set(({ deletions }) => ({
      deletions: deletions.map((candidate) => (candidate === deletion ? { ...candidate, sent: true } : candidate)),
    }));

    const results = await Promise.allSettled(deletion.ids.map(async (id) => remove(id)));
    const removed = deletion.ids.filter((_, index) => results[index]?.status === "fulfilled");

    const failed = results.flatMap((result, index) => {
      const id = deletion.ids[index];

      if (result.status !== "rejected" || id === undefined) {
        return [];
      }

      const error: unknown = result.reason;

      return [{ id, error }];
    });

    set((state) => ({
      deletions: state.deletions.filter((candidate) => candidate.token !== token),
      removed: [...state.removed, ...removed],
    }));

    return { removed, failed };
  },
  forget: (listed) => {
    set(({ removed }) => {
      const kept = removed.filter((id) => listed.has(id));

      return kept.length === removed.length ? {} : { removed: kept };
    });
  },
}));

export function hiddenRunIds({ deletions, removed }: RunDeletionState): ReadonlySet<JobId> {
  return new Set([...deletions.flatMap((deletion) => deletion.ids), ...removed]);
}

export function useHiddenRunIds(): ReadonlySet<JobId> {
  const deletions = useRunDeletionStore((state) => state.deletions);
  const removed = useRunDeletionStore((state) => state.removed);

  return useMemo(() => hiddenRunIds({ deletions, removed }), [deletions, removed]);
}
