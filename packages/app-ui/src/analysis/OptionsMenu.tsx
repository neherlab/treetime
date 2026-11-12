import { type AppCommand, errorMessage, type RunSummary } from "@neherlab/app-contracts";
import { runsSettings } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "@neherlab/design/ui/dialog";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@neherlab/design/ui/dropdown-menu";
import { Kbd } from "@neherlab/design/ui/kbd";
import { toast } from "@neherlab/design/ui/toast";
import { DateTime } from "luxon";
import { useCallback, useMemo, useState } from "react";
import { useDropzone } from "react-dropzone";
import Ellipsis from "~icons/lucide/ellipsis";

import { useApiContext } from "../api/context";
import { headlineText } from "../format";
import { COMMAND_LINE_HOTKEY } from "../hotkeys";
import { withRunSettings } from "../settings/applyRun";
import { commandSettings } from "../settings/catalog";
import { hasEnded } from "../shell/runActions";
import { StatusIcon } from "../shell/StatusIcon";
import { useListedRuns } from "../shell/useListedRuns";
import { useOpenFiles } from "../start/useOpenFiles";
import { useDraftStore } from "../store/draft";
import { dropzoneAccept } from "./fileAccept";

const CONFIG_EXTENSIONS = ["yaml", "yml"];

export function OptionsMenu({ command, onCommandLine }: { command: AppCommand; onCommandLine: () => void }) {
  const [applyOpen, setApplyOpen] = useState(false);
  const openApply = useCallback(() => setApplyOpen(true), []);
  const openFiles = useOpenFiles();

  const onDropAccepted = useCallback((files: File[]) => void openFiles(files), [openFiles]);

  const { getInputProps, open: chooseConfig } = useDropzone({
    noClick: true,
    noDrag: true,
    noKeyboard: true,
    multiple: false,
    onDropAccepted,
    ...dropzoneAccept(CONFIG_EXTENSIONS),
  });

  return (
    <>
      <input {...getInputProps()} aria-label="YAML config to load" />
      <DropdownMenu>
        <DropdownMenuTrigger
          render={<Button type="button" variant="ghost" size="icon-sm" />}
          aria-label="More options actions"
        >
          <Ellipsis aria-hidden />
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" className="w-64">
          <DropdownMenuItem onClick={onCommandLine} className="justify-between">
            Command line
            <Kbd>{COMMAND_LINE_HOTKEY}</Kbd>
          </DropdownMenuItem>
          <DropdownMenuSeparator />
          <DropdownMenuItem onClick={openApply}>Apply settings from a run...</DropdownMenuItem>
          <DropdownMenuItem onClick={chooseConfig}>Load config (YAML)...</DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
      <ApplyRunDialog command={command} open={applyOpen} onOpenChange={setApplyOpen} />
    </>
  );
}

function ApplyRunDialog({
  command,
  open,
  onOpenChange,
}: {
  command: AppCommand;
  open: boolean;
  onOpenChange: (open: boolean) => void;
}) {
  const { client } = useApiContext();
  const { all } = useListedRuns();
  const [busy, setBusy] = useState(false);

  const runs = useMemo(() => all.filter((run) => run.command === command && hasEnded(run.status)), [all, command]);

  const apply = useCallback(
    async (run: RunSummary) => {
      setBusy(true);

      try {
        const { data: settings } = await runsSettings({ client, path: { id: run.id }, throwOnError: true });
        const { draft, load } = useDraftStore.getState();
        const specs = commandSettings(draft.command).settings;

        load({ ...draft, config: withRunSettings(specs, draft.config, run.changed_settings, settings) });
        toast.add({ title: `Applied the settings of ${run.title}` });
        onOpenChange(false);
      } catch (error: unknown) {
        toast.add({ title: "The settings of the run cannot be loaded", description: errorMessage(error) });
      } finally {
        setBusy(false);
      }
    },
    [client, onOpenChange],
  );

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-h-[calc(100%-4rem)] grid-rows-[auto_minmax(0,1fr)] gap-4 sm:max-w-xl">
        <DialogHeader>
          <DialogTitle>Apply settings from a run</DialogTitle>
          <DialogDescription>
            The settings that the run changed from the defaults replace those of the options. The inputs stay.
          </DialogDescription>
        </DialogHeader>
        <div className="min-h-0 overflow-y-auto overscroll-contain">
          {runs.length === 0 ? (
            <p className="text-muted-foreground">No ended run of {commandSettings(command).title}.</p>
          ) : (
            <ul className="grid gap-1">
              {runs.map((run) => (
                <RunChoice key={run.id} run={run} disabled={busy} onPick={apply} />
              ))}
            </ul>
          )}
        </div>
      </DialogContent>
    </Dialog>
  );
}

function RunChoice({
  run,
  disabled,
  onPick,
}: {
  run: RunSummary;
  disabled: boolean;
  onPick: (run: RunSummary) => Promise<void>;
}) {
  const onClick = useCallback(() => void onPick(run), [onPick, run]);
  const created = DateTime.fromISO(run.created_at);
  const changed = run.changed_settings.length;
  const headline = headlineText(run.headline);

  return (
    <li>
      <Button
        type="button"
        variant="ghost"
        className="h-auto w-full justify-start py-1.5 text-left"
        disabled={disabled}
        onClick={onClick}
      >
        <StatusIcon status={run.status} />
        <span className="grid min-w-0 flex-1 gap-0.5">
          <span className="truncate font-bold">{run.title}</span>
          <span className="text-muted-foreground truncate text-xs font-normal">
            {[
              created.isValid ? created.toFormat("d LLL yyyy, HH:mm") : run.created_at,
              `${changed} changed ${changed === 1 ? "setting" : "settings"}`,
              headline,
            ]
              .filter((part) => part !== "")
              .join(", ")}
          </span>
        </span>
      </Button>
    </li>
  );
}
