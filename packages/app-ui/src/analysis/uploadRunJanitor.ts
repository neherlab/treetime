import { errorMessage, type JobId } from "@neherlab/app-contracts";
import { ApiError, type ApiClient, runsDelete, runsGet } from "@neherlab/app-contracts/client";

import {
  heldCreations,
  type HeldCreations,
  heldUploadRuns,
  type HeldDrafts,
  releasedUploadRuns,
} from "../store/draftRules";

export class UploadRunJanitor {
  readonly #store: JanitorStore;
  readonly #retire: (runId: JobId) => Promise<void>;
  readonly #report: (problem: JanitorProblem) => void;
  readonly #retiring = new Set<JobId>();
  readonly #pending = new Set<Promise<void>>();
  readonly #unsubscribe: () => void;

  constructor(store: JanitorStore, retire: (runId: JobId) => Promise<void>, report: (problem: JanitorProblem) => void) {
    this.#store = store;
    this.#retire = retire;
    this.#report = report;
    this.#unsubscribe = store.subscribe((after, before) => {
      for (const runId of releasedUploadRuns(before, after)) {
        this.#track(this.#retireWhenIdle(runId));
      }

      const kept = heldCreations(after);

      for (const creation of heldCreations(before)) {
        if (!kept.has(creation)) {
          this.#track(this.#retireCreated(creation));
        }
      }
    });
  }

  dispose(): void {
    this.#unsubscribe();
  }

  async settled(): Promise<void> {
    while (this.#pending.size > 0) {
      await Promise.all(this.#pending);
    }
  }

  #track(work: Promise<void>): void {
    this.#pending.add(work);
    void work.finally(() => {
      this.#pending.delete(work);
    });
  }

  async #retireCreated(creation: Promise<JobId>): Promise<void> {
    try {
      await this.#retireWhenIdle(await creation);
    } catch (error: unknown) {
      this.#report({
        title: "The upload run of a discarded form was not created",
        description: errorMessage(error),
      });
    }
  }

  async #retireWhenIdle(runId: JobId): Promise<void> {
    if (this.#retiring.has(runId)) {
      return;
    }

    this.#retiring.add(runId);

    try {
      await this.#idle(runId);

      if (!heldUploadRuns(this.#store.getState()).has(runId)) {
        await this.#retire(runId);
      }
    } catch (error: unknown) {
      this.#report({ title: "The uploads of a discarded form remain", description: errorMessage(error) });
    } finally {
      this.#retiring.delete(runId);
    }
  }

  async #idle(runId: JobId): Promise<void> {
    const { promise, resolve } = Promise.withResolvers<undefined>();

    const check = (state: JanitorState) => {
      if (!state.busyRuns.includes(runId)) {
        unsubscribe();
        resolve(undefined);
      }
    };

    const unsubscribe = this.#store.subscribe(check);
    check(this.#store.getState());

    return promise;
  }
}

export async function deleteUnstartedRun(client: ApiClient, runId: JobId): Promise<void> {
  try {
    const { data: record } = await runsGet({ client, path: { id: runId }, throwOnError: true });

    if (record.status === "created") {
      await runsDelete({ client, path: { id: runId }, throwOnError: true });
    }
  } catch (error: unknown) {
    if (!(error instanceof ApiError && error.response.code === "not_found")) {
      throw error;
    }
  }
}

export interface JanitorProblem {
  title: string;
  description: string;
}

export interface JanitorState extends HeldDrafts, HeldCreations {
  busyRuns: readonly JobId[];
}

export interface JanitorStore {
  getState: () => JanitorState;
  subscribe: (listener: (state: JanitorState, previous: JanitorState) => void) => () => void;
}
