import type { JobId, RunRecord, UiDraft } from "@neherlab/app-contracts";
import { type ApiClient, storedInputs } from "@neherlab/app-contracts/client";

import { baseName } from "../settings/inputs";
import {
  draftInputPaths,
  forgetRuns,
  type ForgottenRuns,
  type LostInputs,
  type QueuedExample,
  referencedRuns,
} from "../store/draftSchema";
import { findRun } from "./pendingUpload";

export async function forgetLostRuns(client: ApiClient, draft: UiDraft): Promise<CheckedRuns> {
  const { lost, runs } = await lostInputs(client, [draft]);

  return { ...forgetRuns(draft, lost), runs };
}

export async function lostInputs(client: ApiClient, drafts: readonly UiDraft[]): Promise<CheckedInputs> {
  const ids = [...new Set(drafts.flatMap(referencedRuns))];
  const paths = [...new Set(drafts.flatMap(draftInputPaths))];
  const records = await Promise.all(ids.map(async (id) => findRun(client, id)));
  const runs = new Map(records.flatMap((record) => (record === undefined ? [] : [[record.id, record] as const])));

  const { data: stored } =
    paths.length === 0
      ? { data: { missing: [] } }
      : await storedInputs({ client, body: { paths }, throwOnError: true });

  return { lost: { runs: new Set(ids.filter((id) => !runs.has(id))), paths: new Set(stored.missing) }, runs };
}

export function reloadedFilesMessage(labels: readonly string[]) {
  return {
    title: "The files of the form belonged to the previous page; choose them again",
    description: labels.join(", "),
  };
}

export function refetchedExamplesMessage(examples: readonly QueuedExample[]) {
  return {
    title: "The example files are being fetched again",
    description: `The back end no longer has ${examples.map((example) => baseName(example.path)).join(", ")}. Run again when they are back in the form.`,
  };
}

export function lostFilesMessage(labels: readonly string[]) {
  return {
    title: "Some input files are no longer available",
    description: `The back end no longer has ${labels.join(", ")}, so the form left them out. Add them again.`,
  };
}

export interface CheckedRuns extends ForgottenRuns {
  runs: ReadonlyMap<JobId, RunRecord>;
}

export interface CheckedInputs {
  lost: LostInputs;
  runs: ReadonlyMap<JobId, RunRecord>;
}
