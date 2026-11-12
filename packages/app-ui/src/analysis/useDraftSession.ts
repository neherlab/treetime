import { toast } from "@neherlab/design/ui/toast";
import { useNavigate, useParams } from "@tanstack/react-router";
import { useCallback, useEffect, useRef } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import type { DiscardReason } from "../store/draftRules";
import { deleteUnstartedRun, UploadRunJanitor } from "./uploadRunJanitor";
import { useQueuedExamples } from "./useExampleData";

const UNDO_TIMEOUT_MS = 10_000;

const UNDO_TITLES: Readonly<Record<DiscardReason, string>> = {
  reset: "The form was reset",
  new: "A new analysis replaced your own draft",
  replaced: "Opening the run set the edited variant aside",
  revert: "The settings of the run replaced the edited options",
};

export function useDraftSession(): void {
  useUploadRunJanitor();
  useDiscardUndo();
  useQueuedExamples();
}

export function useNewAnalysis(): () => void {
  const navigate = useNavigate();
  const discard = useDraftStore((state) => state.discard);

  return useCallback(() => {
    discard("new", "timetree");
    void navigate({ to: "/new" });
  }, [discard, navigate]);
}

function useUploadRunJanitor(): void {
  const { client } = useApiContext();

  useEffect(() => {
    const janitor = new UploadRunJanitor(
      useDraftStore,
      async (runId) => deleteUnstartedRun(client, runId),
      (problem) => {
        toast.add(problem);
      },
    );

    return () => {
      janitor.dispose();
    };
  }, [client]);
}

function useDiscardUndo(): void {
  const undo = useDraftStore((state) => state.undo);
  const { id: shownRunId } = useParams({ strict: false });
  const shown = useRef(shownRunId);

  useEffect(() => {
    shown.current = shownRunId;
  }, [shownRunId]);

  useEffect(() => {
    if (undo === undefined) {
      return undefined;
    }

    const { token, reason } = undo;

    const id = toast.add({
      title: UNDO_TITLES[reason],
      timeout: UNDO_TIMEOUT_MS,
      actionProps: {
        children: "Undo",
        onClick: () => {
          useDraftStore.getState().undoDiscard(token, shown.current);
          toast.close(id);
        },
      },
      onClose: () => {
        useDraftStore.getState().dropUndo(token);
      },
    });

    return () => {
      toast.close(id);
    };
  }, [undo]);
}
