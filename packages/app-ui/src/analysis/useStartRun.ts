import { type SparseConfig, errorMessage, type AppCommand, type StartRunRequest } from "@neherlab/app-contracts";
import { runsCreate, runsStart } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import { exampleDatasetName } from "./datasetLabel";
import { forgetLostRuns, lostFilesMessage } from "./lostRuns";
import { pendingUpload } from "./pendingUpload";
import { whileRunBusy } from "./runWork";

export function useStartRun(command: AppCommand) {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const toasts = useToastManager();

  return useCallback(
    async (config: SparseConfig) => {
      const { draft, finishStart, load, queueExamples } = useDraftStore.getState();

      try {
        const kept = await forgetLostRuns(client, draft);

        if (kept.dropped.length > 0 || kept.examples.length > 0) {
          load(kept.draft);
          queueExamples(kept.examples);

          if (kept.dropped.length > 0) {
            toasts.add(lostFilesMessage(kept.dropped));
          }

          return;
        }

        const uploadRunId = kept.draft.upload_run_id;
        const upload = pendingUpload(uploadRunId === undefined ? undefined : kept.runs.get(uploadRunId));

        const example = exampleDatasetName(kept.draft);
        const source = kept.draft.from_run_id;
        const named: Pick<StartRunRequest, "example" | "source_run_id"> = {};

        if (example !== undefined) {
          named.example = example;
        }

        if (source !== undefined) {
          named.source_run_id = source;
        }

        const { data: record } =
          upload === undefined
            ? await runsCreate({ client, body: { command, config, defer_start: false, ...named }, throwOnError: true })
            : await whileRunBusy(upload.id, async () =>
                runsStart({ client, path: { id: upload.id }, body: { command, config, ...named }, throwOnError: true }),
              );

        finishStart();
        await navigate({ to: "/runs/$id/results", params: { id: record.id } });
      } catch (error: unknown) {
        toasts.add({
          title: "The run cannot be started",
          description: errorMessage(error),
        });
      }
    },
    [client, command, navigate, toasts],
  );
}
