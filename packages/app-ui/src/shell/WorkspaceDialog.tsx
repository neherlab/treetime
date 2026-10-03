import { errorMessage, type Workspace, type WorkspaceShell } from "@neherlab/app-contracts";
import { runsList, workspace as getWorkspace, workspaceUpdate } from "@neherlab/app-contracts/client";
import { useQueryClient } from "@tanstack/react-query";
import { useCallback, useState } from "react";
import FolderOpen from "~icons/lucide/folder-open";

import { useApiContext } from "../api/context";
import { useApi } from "../api/hooks";
import { useWorkspaceShell } from "../platform";
import { useShellStore } from "../store/shell";
import { Alert, AlertDescription } from "../ui/alert";
import { Button } from "../ui/button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "../ui/dialog";
import { useToastManager } from "../ui/toast";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const RUNS_FOLDER = "Runs folder";

export function WorkspaceButton() {
  const shell = useWorkspaceShell();
  const setOpen = useShellStore((state) => state.setWorkspaceOpen);
  const openDialog = useCallback(() => setOpen(true), [setOpen]);

  if (shell === null) {
    return null;
  }

  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button variant="ghost" size="icon" aria-label={RUNS_FOLDER} onClick={openDialog}>
            <FolderOpen aria-hidden />
          </Button>
        }
      />
      <TooltipContent>{RUNS_FOLDER}: where TreeTime keeps the runs. Click to change.</TooltipContent>
    </Tooltip>
  );
}

export function WorkspaceDialog() {
  const shell = useWorkspaceShell();
  const open = useShellStore((state) => state.workspaceOpen);
  const setOpen = useShellStore((state) => state.setWorkspaceOpen);
  const close = useCallback(() => setOpen(false), [setOpen]);

  if (shell === null) {
    return null;
  }

  return (
    <Dialog open={open} onOpenChange={setOpen}>
      <DialogContent>{open && <WorkspaceForm shell={shell} onChanged={close} />}</DialogContent>
    </Dialog>
  );
}

function WorkspaceForm({ shell, onChanged }: { shell: WorkspaceShell; onChanged: () => void }) {
  const { data: workspace } = useApi((context) => getWorkspace(context));
  const { data: runList } = useApi((context) => runsList(context));
  const change = useWorkspaceChange(shell, onChanged);
  const [busy, setBusy] = useState(false);

  const apply = useCallback(
    async (path: string | null) => {
      setBusy(true);
      await change(path);
      setBusy(false);
    },
    [change],
  );

  const choose = useCallback(async () => {
    const path = await shell.pickFolder({ title: RUNS_FOLDER });

    if (path !== null) {
      await apply(path);
    }
  }, [apply, shell]);

  const onChoose = useCallback(() => void choose(), [choose]);
  const onResetToDefault = useCallback(() => void apply(null), [apply]);
  const computing = runList?.active_runs ?? 0;

  return (
    <>
      <DialogHeader>
        <DialogTitle>{RUNS_FOLDER}</DialogTitle>
        <DialogDescription>
          TreeTime keeps each run, with its inputs and results, in this folder. Runs in another folder appear when you
          switch back to it.
        </DialogDescription>
      </DialogHeader>
      <WorkspacePath workspace={workspace} />
      {computing > 0 && (
        <Alert>
          <AlertDescription>
            {computing === 1 ? "1 run is computing" : `${computing} runs are computing`}. Changing the folder stops
            them, and they stay in the current folder as interrupted runs.
          </AlertDescription>
        </Alert>
      )}
      <div className="flex flex-wrap justify-end gap-2">
        <Button
          variant="outline"
          disabled={busy || workspace === undefined || workspace.path === workspace.default_path}
          onClick={onResetToDefault}
        >
          Use the default folder
        </Button>
        <Button disabled={busy} onClick={onChoose}>
          <FolderOpen aria-hidden />
          Choose a folder
        </Button>
      </div>
    </>
  );
}

function WorkspacePath({ workspace }: { workspace: Workspace | undefined }) {
  if (workspace === undefined) {
    return null;
  }

  return (
    <dl className="grid gap-1 text-sm">
      <dt className="text-muted-foreground">Current folder</dt>
      <dd className="font-mono break-all">{workspace.path}</dd>
      {workspace.path !== workspace.default_path && (
        <>
          <dt className="text-muted-foreground">Default folder</dt>
          <dd className="font-mono break-all">{workspace.default_path}</dd>
        </>
      )}
    </dl>
  );
}

function useWorkspaceChange(shell: WorkspaceShell, onChanged: () => void) {
  const { client } = useApiContext();
  const queryClient = useQueryClient();
  const toasts = useToastManager();

  return useCallback(
    async (path: string | null) => {
      try {
        await workspaceUpdate({ client, body: { path }, throwOnError: true });
        await shell.restartBackend();
        await queryClient.resetQueries();
        onChanged();
        toasts.add({ title: "The runs folder changed", description: path ?? "TreeTime uses its default folder." });
      } catch (error: unknown) {
        toasts.add({ title: "The runs folder cannot be changed", description: errorMessage(error) });
      }
    },
    [client, onChanged, queryClient, shell, toasts],
  );
}
