import type { RunRecord } from "@neherlab/app-contracts";
import { runsFiles } from "@neherlab/app-contracts/client";
import { Alert, AlertDescription } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "@neherlab/design/ui/dialog";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@neherlab/design/ui/dropdown-menu";
import { useCallback, useMemo, useState } from "react";
import ChevronDown from "~icons/lucide/chevron-down";
import Download from "~icons/lucide/download";

import { useApi } from "../api/hooks";
import { formatBytes } from "../format";
import { type DownloadEntry, downloadEntries, totalSize } from "../results/files";
import { OutputFilesTable, useSaveRunOutput } from "./OutputFiles";

export function DownloadMenu({ record, packageName }: { record: RunRecord; packageName: string | undefined }) {
  const save = useSaveRunOutput();
  const [allOpen, setAllOpen] = useState(false);

  const { data: files, error } = useApi((context) => runsFiles({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const entries = useMemo(() => (files === undefined ? undefined : downloadEntries(files)), [files]);

  const onEntry = useCallback(
    (entry: DownloadEntry) => {
      switch (entry.kind) {
        case "package":
          if (packageName !== undefined) {
            void save({ id: record.id, name: packageName }, "The run package cannot be downloaded");
          }

          return;
        case "file":
          void save(
            { id: record.id, path: entry.file.path, name: entry.file.download_name },
            "The file cannot be downloaded",
          );

          return;
        case "all":
          setAllOpen(true);
      }
    },
    [packageName, record.id, save],
  );

  return (
    <>
      <DropdownMenu>
        <DropdownMenuTrigger render={<Button type="button" variant="outline" size="sm" />}>
          <Download aria-hidden />
          Download
          <ChevronDown aria-hidden />
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" className="w-80">
          {error !== null && (
            <DropdownMenuItem disabled>The file list cannot be loaded: {error.message}</DropdownMenuItem>
          )}
          {entries === undefined && error === null && <DropdownMenuItem disabled>Loading the files</DropdownMenuItem>}
          {entries?.map((entry) => (
            <DownloadItem
              key={entryKey(entry)}
              entry={entry}
              disabled={entry.kind === "package" && packageName === undefined}
              onSelect={onEntry}
            />
          ))}
        </DropdownMenuContent>
      </DropdownMenu>
      <Dialog open={allOpen} onOpenChange={setAllOpen}>
        <DialogContent className="max-h-[calc(100%-4rem)] grid-rows-[auto_minmax(0,1fr)] gap-4 sm:max-w-4xl">
          <DialogHeader>
            <DialogTitle>All files</DialogTitle>
            <DialogDescription>
              {files === undefined
                ? record.title
                : `${files.length} files of ${record.title}, ${formatBytes(totalSize(files))}`}
            </DialogDescription>
          </DialogHeader>
          <div className="min-h-0 overflow-y-auto overscroll-contain">
            {error !== null && (
              <Alert variant="destructive">
                <AlertDescription>The file list cannot be loaded: {error.message}</AlertDescription>
              </Alert>
            )}
            {files !== undefined && <OutputFilesTable runId={record.id} files={files} />}
          </div>
        </DialogContent>
      </Dialog>
    </>
  );
}

function DownloadItem({
  entry,
  disabled,
  onSelect,
}: {
  entry: DownloadEntry;
  disabled: boolean;
  onSelect: (entry: DownloadEntry) => void;
}) {
  const onClick = useCallback(() => onSelect(entry), [entry, onSelect]);

  if (entry.kind === "package") {
    return (
      <>
        <DropdownMenuItem onClick={onClick} disabled={disabled} className="grid gap-0.5">
          <span className="font-bold">Run package (.zip)</span>
          <span className="text-muted-foreground text-xs">Inputs, settings, log, and outputs</span>
        </DropdownMenuItem>
        <DropdownMenuSeparator />
      </>
    );
  }

  if (entry.kind === "file") {
    return (
      <DropdownMenuItem onClick={onClick} className="justify-between gap-3">
        <span>{entry.file.label}</span>
        <span className="text-muted-foreground text-xs tabular-nums">{formatBytes(entry.file.size)}</span>
      </DropdownMenuItem>
    );
  }

  return (
    <>
      <DropdownMenuSeparator />
      <DropdownMenuItem onClick={onClick}>All files...</DropdownMenuItem>
    </>
  );
}

function entryKey(entry: DownloadEntry): string {
  return entry.kind === "file" ? `file:${entry.file.path}` : entry.kind;
}
