import { errorMessage, type JobId, type RunStatus, type RunSummary } from "@neherlab/app-contracts";
import { type ApiClient, runsList, runsUpdate } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "@neherlab/design/ui/dialog";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@neherlab/design/ui/dropdown-menu";
import { Input } from "@neherlab/design/ui/input";
import { toast } from "@neherlab/design/ui/toast";
import { useNavigate, useParams } from "@tanstack/react-router";
import { useCallback, useMemo, useState } from "react";
import Ellipsis from "~icons/lucide/ellipsis";

import { lostFilesMessage, lostInputs } from "../analysis/lostRuns";
import { useApiContext } from "../api/context";
import { useApi, useApiMutation } from "../api/hooks";
import { useSaveRunOutput } from "../runs/OutputFiles";
import { baseName } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { ownDraft } from "../store/draftRules";
import { forgetRuns, runSourcePaths } from "../store/draftSchema";
import { useRunDeletionStore } from "../store/runDeletion";
import { useShellStore } from "../store/shell";
import { onPageHide } from "./pageHide";
import { removeRun } from "./removeRun";

const UNDO_TIMEOUT_MS = 10_000;

const CANCEL_WAIT_MS = 60_000;

const ENDED: Readonly<Record<RunStatus, boolean>> = {
  created: false,
  running: false,
  ok: true,
  error: true,
  cancelled: true,
  interrupted: true,
};

export type ActionRun = Pick<RunSummary, "id" | "title" | "pinned" | "status" | "source_run_id" | "package_name">;

interface RunActions {
  rename: (run: ActionRun) => void;
  setPinned: (run: ActionRun, pinned: boolean) => Promise<void>;
  downloadPackage: (run: ActionRun) => Promise<void>;
  remove: (runs: readonly ActionRun[]) => void;
  compareWithSource: (run: ActionRun) => void;
}

export function hasEnded(status: RunStatus): boolean {
  return ENDED[status];
}

export function useRunActions(): RunActions {
  const navigate = useNavigate();
  const save = useSaveRunOutput();
  const deleteRuns = useDeleteRuns();
  const setRenamingId = useShellStore((state) => state.setRenamingId);
  const setConfirmingDelete = useShellStore((state) => state.setConfirmingDelete);

  const { mutateAsync: update } = useApiMutation((context, { id, pinned }: { id: JobId; pinned: boolean }) =>
    runsUpdate({ ...context, path: { id }, body: { pinned } }),
  );

  const setPinned = useCallback(
    async (run: ActionRun, pinned: boolean) => {
      try {
        await update({ id: run.id, pinned });
      } catch (error: unknown) {
        toast.add({
          title: `${run.title} cannot be ${pinned ? "pinned" : "unpinned"}`,
          description: errorMessage(error),
        });
      }
    },
    [update],
  );

  const downloadPackage = useCallback(
    async (run: ActionRun) => save({ id: run.id, name: run.package_name }, "The run package cannot be downloaded"),
    [save],
  );

  const remove = useCallback(
    (runs: readonly ActionRun[]) => {
      const ids = runs.map((run) => run.id);

      if (runs.some((run) => run.status === "running") || draftFilesOf(ids).length > 0) {
        setConfirmingDelete(ids);
      } else {
        deleteRuns(ids);
      }
    },
    [deleteRuns, setConfirmingDelete],
  );

  const compareWithSource = useCallback(
    (run: ActionRun) => {
      if (run.source_run_id !== undefined) {
        void navigate({ to: "/compare/$a/$b", params: { a: run.source_run_id, b: run.id } });
      }
    },
    [navigate],
  );

  const rename = useCallback((run: ActionRun) => setRenamingId(run.id), [setRenamingId]);

  return useMemo(
    () => ({ rename, setPinned, downloadPackage, remove, compareWithSource }),
    [compareWithSource, downloadPackage, remove, rename, setPinned],
  );
}

export function RunActionsMenu({ run, children }: { run: ActionRun; children?: React.ReactNode }) {
  const actions = useRunActions();
  const { data } = useApi((context) => runsList(context));
  const source = data?.runs.find((candidate) => candidate.id === run.source_run_id);
  const ended = ENDED[run.status];

  const onRename = useCallback(() => actions.rename(run), [actions, run]);
  const onPin = useCallback(() => void actions.setPinned(run, !run.pinned), [actions, run]);
  const onCompare = useCallback(() => actions.compareWithSource(run), [actions, run]);
  const onDownload = useCallback(() => void actions.downloadPackage(run), [actions, run]);
  const onDelete = useCallback(() => actions.remove([run]), [actions, run]);

  return (
    <DropdownMenu>
      <DropdownMenuTrigger
        render={<Button type="button" variant="ghost" size="icon-xs" />}
        aria-label={`Actions of ${run.title}`}
      >
        <Ellipsis aria-hidden />
      </DropdownMenuTrigger>
      <DropdownMenuContent align="end" className="w-52">
        {children}
        <DropdownMenuItem onClick={onRename}>Rename</DropdownMenuItem>
        <DropdownMenuItem onClick={onPin}>{run.pinned ? "Unpin" : "Pin"}</DropdownMenuItem>
        <DropdownMenuItem
          onClick={onCompare}
          disabled={run.status !== "ok" || source === undefined || source.status !== "ok"}
        >
          Compare with source
        </DropdownMenuItem>
        <DropdownMenuItem onClick={onDownload} disabled={!ended}>
          Download run package
        </DropdownMenuItem>
        <DropdownMenuItem onClick={onDelete} className="text-destructive focus:text-destructive">
          Delete
        </DropdownMenuItem>
      </DropdownMenuContent>
    </DropdownMenu>
  );
}

export function RunActionDialogs() {
  return (
    <>
      <RenameRunDialog />
      <ConfirmDeleteDialog />
    </>
  );
}

function RenameRunDialog() {
  const renamingId = useShellStore((state) => state.renamingId);
  const setRenamingId = useShellStore((state) => state.setRenamingId);
  const { data } = useApi((context) => runsList(context));
  const run = data?.runs.find((candidate) => candidate.id === renamingId);
  const close = useCallback(() => setRenamingId(undefined), [setRenamingId]);

  const onOpenChange = useCallback(
    (open: boolean) => {
      if (!open) {
        close();
      }
    },
    [close],
  );

  return (
    <Dialog open={run !== undefined} onOpenChange={onOpenChange}>
      <DialogContent>
        <DialogHeader>
          <DialogTitle>Rename run</DialogTitle>
        </DialogHeader>
        {run !== undefined && <RenameForm key={run.id} run={run} close={close} />}
      </DialogContent>
    </Dialog>
  );
}

function RenameForm({ run, close }: { run: ActionRun; close: () => void }) {
  const [title, setTitle] = useState(run.title);

  const { mutateAsync: rename, isPending } = useApiMutation((context, next: string) =>
    runsUpdate({ ...context, path: { id: run.id }, body: { title: next } }),
  );

  const onChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => setTitle(event.target.value), []);

  const submit = useCallback(async () => {
    const next = title.trim();

    if (next !== "" && next !== run.title) {
      try {
        await rename(next);
      } catch (error: unknown) {
        toast.add({ title: "The run cannot be renamed", description: errorMessage(error) });

        return;
      }
    }

    close();
  }, [close, rename, run.title, title]);

  const onSubmit = useCallback(
    (event: React.SubmitEvent) => {
      event.preventDefault();
      void submit();
    },
    [submit],
  );

  const focusInput = useCallback((input: HTMLInputElement | null) => {
    input?.focus();
    input?.select();
  }, []);

  return (
    <form className="grid gap-4" onSubmit={onSubmit}>
      <Input type="text" aria-label="Name of the run" value={title} onChange={onChange} ref={focusInput} />
      <div className="flex justify-end gap-2">
        <Button type="button" variant="ghost" onClick={close}>
          Cancel
        </Button>
        <Button type="submit" disabled={isPending || title.trim() === ""}>
          Rename
        </Button>
      </div>
    </form>
  );
}

function ConfirmDeleteDialog() {
  const ids = useShellStore((state) => state.confirmingDelete);
  const setConfirmingDelete = useShellStore((state) => state.setConfirmingDelete);
  const deleteRuns = useDeleteRuns();
  const { data } = useApi((context) => runsList(context));
  const close = useCallback(() => setConfirmingDelete(undefined), [setConfirmingDelete]);

  const onOpenChange = useCallback(
    (open: boolean) => {
      if (!open) {
        close();
      }
    },
    [close],
  );

  const confirm = useCallback(() => {
    if (ids !== undefined) {
      deleteRuns(ids);
    }

    close();
  }, [close, deleteRuns, ids]);

  const runs = (data?.runs ?? []).filter((run) => ids?.includes(run.id) === true);
  const computing = runs.filter((run) => run.status === "running");
  const draftFiles = ids === undefined ? [] : draftFilesOf(ids);

  return (
    <Dialog open={ids !== undefined} onOpenChange={onOpenChange}>
      <DialogContent>
        <DialogHeader>
          <DialogTitle>
            {runs.length === 1 ? `Delete ${runs[0]?.title ?? "the run"}?` : `Delete ${runs.length} runs?`}
          </DialogTitle>
          {computing.length > 0 && (
            <DialogDescription>
              {computing.length === 1 ? "The run is" : `${computing.length} runs are`} still computing. Deleting cancels
              the computation first.
            </DialogDescription>
          )}
          {draftFiles.length > 0 && (
            <DialogDescription>
              The open analysis reads {draftFiles.join(", ")} from {runs.length === 1 ? "this run" : "these runs"}. The
              form leaves these files out after the delete; add them again to run the analysis.
            </DialogDescription>
          )}
        </DialogHeader>
        <div className="flex justify-end gap-2">
          <Button type="button" variant="ghost" onClick={close}>
            Cancel
          </Button>
          <Button type="button" variant="destructive" onClick={confirm}>
            Delete
          </Button>
        </div>
      </DialogContent>
    </Dialog>
  );
}

function useDeleteRuns(): (ids: readonly JobId[]) => void {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const shownId = useParams({ strict: false, select: (params) => params.id });
  const { data: list } = useApi((context) => runsList(context));

  return useCallback(
    (ids: readonly JobId[]) => {
      const { hide, undo } = useRunDeletionStore.getState();
      const token = hide(ids);

      if (token === undefined) {
        return;
      }

      const { selectedIds, clearSelection } = useShellStore.getState();

      if (selectedIds.some((id) => ids.includes(id))) {
        clearSelection();
      }

      if (shownId !== undefined && ids.includes(shownId)) {
        void navigate({ to: "/" });
      }

      const titles = new Map((list?.runs ?? []).map((run) => [run.id, run.title]));

      const commit = () => {
        stopOnPageHide();
        void commitDeletion(client, token, titles, async (id) =>
          removeRun(client, id, AbortSignal.timeout(CANCEL_WAIT_MS)),
        );
      };

      const stopOnPageHide = onPageHide(window, commit);

      const toastId = toast.add({
        title: ids.length === 1 ? "The run was deleted" : `${ids.length} runs were deleted`,
        timeout: UNDO_TIMEOUT_MS,
        actionProps: {
          children: "Undo",
          onClick: () => {
            stopOnPageHide();
            undo(token);
            toast.close(toastId);
          },
        },
        onClose: commit,
      });
    },
    [client, list, navigate, shownId],
  );
}

async function commitDeletion(
  client: ApiClient,
  token: symbol,
  titles: ReadonlyMap<JobId, string>,
  remove: (id: JobId) => Promise<void>,
): Promise<void> {
  const { removed, failed } = await useRunDeletionStore.getState().commit(token, remove);

  for (const { id, error } of failed) {
    toast.add({ title: `${titles.get(id) ?? "The run"} cannot be deleted`, description: errorMessage(error) });
  }

  if (removed.length === 0) {
    return;
  }

  try {
    await forgetDeletedRuns(client);
  } catch (error: unknown) {
    toast.add({ title: "The files of the form cannot be checked", description: errorMessage(error) });
  }
}

function draftFilesOf(ids: readonly JobId[]): string[] {
  const draft = ownDraft(useDraftStore.getState());
  const runs = new Set(ids);

  return forgetRuns(draft, { runs, paths: new Set(runSourcePaths(draft, runs)) }).dropped;
}

async function forgetDeletedRuns(client: ApiClient): Promise<void> {
  const { draft, parked } = useDraftStore.getState();
  const drafts = parked === undefined ? [draft] : [draft, parked.draft];
  const { lost } = await lostInputs(client, drafts);
  const forgotten = drafts.map((each) => forgetRuns(each, lost));
  const [current] = forgotten;

  const dropped = [
    ...new Set([
      ...forgotten.flatMap((each) => each.dropped),
      ...forgotten.slice(1).flatMap((each) => each.examples.map((example) => baseName(example.path))),
    ]),
  ];

  useDraftStore.getState().forgetInputs(lost);

  if (current !== undefined && current.examples.length > 0) {
    useDraftStore.getState().queueExamples(current.examples);
  }

  if (dropped.length > 0) {
    toast.add(lostFilesMessage(dropped));
  }
}
