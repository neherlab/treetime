import { errorMessage, type Workspace } from "@neherlab/app-contracts";
import { runsList, workspace as getWorkspace, workspaceUpdate } from "@neherlab/app-contracts/client";
import { useQueryClient } from "@tanstack/react-query";
import { useCallback, useEffect, useRef, useState } from "react";
import FolderOpen from "~icons/lucide/folder-open";

import { useApiContext } from "../api/context";
import { resetApiQueries, useApi } from "../api/hooks";
import type { Host } from "../host";
import { useHost } from "../host-context";
import { useShellStore } from "../store/shell";
import { Alert, AlertDescription } from "../ui/alert";
import { Button } from "../ui/button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "../ui/dialog";
import { useToastManager } from "../ui/toast";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const RUNS_FOLDER = "Runs folder";

export function WorkspaceButton() {
  const host = useHost();
  const setOpen = useShellStore((state) => state.setWorkspaceOpen);
  const openDialog = useCallback(() => setOpen(true), [setOpen]);

  if (host === null) {
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
  const host = useHost();
  const open = useShellStore((state) => state.workspaceOpen);
  const setOpen = useShellStore((state) => state.setWorkspaceOpen);
  const close = useCallback(() => setOpen(false), [setOpen]);

  if (host === null) {
    return null;
  }

  return (
    <Dialog open={open} onOpenChange={setOpen}>
      <WorkspaceErrorNotice />
      <DialogContent>{open && <WorkspaceForm host={host} onChanged={close} />}</DialogContent>
    </Dialog>
  );
}

function WorkspaceErrorNotice() {
  const { data: workspace } = useApi((context) => getWorkspace(context));
  const setOpen = useShellStore((state) => state.setWorkspaceOpen);
  const toasts = useToastManager();
  const shown = useRef(false);
  const error = workspace?.error ?? undefined;

  useEffect(() => {
    if (error === undefined || shown.current) {
      return;
    }

    shown.current = true;
    toasts.add({
      type: "warning",
      title: "The runs folder cannot be opened",
      description: error,
      timeout: 0,
      actionProps: { children: RUNS_FOLDER, onClick: () => setOpen(true) },
    });
  }, [error, setOpen, toasts]);

  return null;
}

function WorkspaceForm({ host, onChanged }: { host: Host; onChanged: () => void }) {
  const { data: workspace } = useApi((context) => getWorkspace(context));
  const { data: runList } = useApi((context) => runsList(context));
  const change = useWorkspaceChange(host, onChanged);
  const [busy, setBusy] = useState(false);

  const apply = useCallback(
    async (path: string | undefined) => {
      setBusy(true);
      await change(path);
      setBusy(false);
    },
    [change],
  );

  const choose = useCallback(async () => {
    const path = await host.pickFolder({ title: RUNS_FOLDER });

    if (path !== undefined) {
      await apply(path);
    }
  }, [apply, host]);

  const onChoose = useCallback(() => void choose(), [choose]);
  const onResetToDefault = useCallback(() => void apply(undefined), [apply]);
  const computing = runList?.active_runs ?? 0;
  const fixedBy = workspace?.fixed_by ?? undefined;
  const locked = busy || fixedBy !== undefined;

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
      {workspace?.error !== undefined && (
        <Alert variant="destructive">
          <AlertDescription>{workspace.error}</AlertDescription>
        </Alert>
      )}
      {fixedBy !== undefined && (
        <Alert>
          <AlertDescription>
            The environment variable <code className="font-mono">{fixedBy}</code> sets this folder. Unset it and restart
            TreeTime to choose the folder here.
          </AlertDescription>
        </Alert>
      )}
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
          disabled={locked || workspace === undefined || workspace.path === workspace.default_path}
          onClick={onResetToDefault}
        >
          Use the default folder
        </Button>
        <Button disabled={locked} onClick={onChoose}>
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

function useWorkspaceChange(host: Host, onChanged: () => void) {
  const { client } = useApiContext();
  const queryClient = useQueryClient();
  const toasts = useToastManager();

  return useCallback(
    async (path: string | undefined) => {
      try {
        await workspaceUpdate({ client, body: path === undefined ? {} : { path }, throwOnError: true });
        await host.restartBackend();
        await resetApiQueries(queryClient);
        onChanged();
        toasts.add({ title: "The runs folder changed", description: path ?? "TreeTime uses its default folder." });
      } catch (error: unknown) {
        toasts.add({ title: "The runs folder cannot be changed", description: errorMessage(error) });
      }
    },
    [client, onChanged, queryClient, host, toasts],
  );
}
