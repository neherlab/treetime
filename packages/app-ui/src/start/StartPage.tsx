import { errorMessage } from "@neherlab/app-contracts";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { Spinner } from "@neherlab/design/ui/spinner";
import { Link } from "@tanstack/react-router";
import { cn } from "cn";
import { useCallback, useMemo, useState } from "react";
import { useDropzone } from "react-dropzone";
import Upload from "~icons/lucide/upload";

import { headlineText } from "../format";
import { useHost } from "../host-context";
import { withListKeys } from "../results/listKeys";
import { commandSettings } from "../settings/catalog";
import { groupRuns, newestRun, type RunGroup } from "../shell/runList";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { useListedRuns } from "../shell/useListedRuns";
import { ExampleGallery } from "./ExampleGallery";
import { type UnplacedName, useOpenFiles } from "./useOpenFiles";

const RECENT_GROUPS = 4;

const START_COMMAND = "timetree";

export function StartPage() {
  const host = useHost();

  return (
    <div className="mx-auto grid w-full max-w-5xl gap-10 px-6 py-8">
      <section aria-labelledby="start-title" className="grid gap-3">
        <h1 id="start-title" className="text-xl font-bold">
          Start an analysis
        </h1>
        <DropZone />
        {host === null && (
          <p className="text-muted-foreground text-sm">Your files stay in this browser. Nothing is uploaded.</p>
        )}
      </section>
      <RecentWork />
      <ExampleGallery />
    </div>
  );
}

function DropZone() {
  const openFiles = useOpenFiles();
  const [busy, setBusy] = useState(false);
  const [unplaced, setUnplaced] = useState<readonly UnplacedName[]>([]);
  const [failure, setFailure] = useState<string | undefined>(undefined);

  const open = useCallback(
    async (files: readonly File[]) => {
      setBusy(true);
      setFailure(undefined);

      try {
        setUnplaced((await openFiles(files)).unplaced);
      } catch (error: unknown) {
        setUnplaced([]);
        setFailure(errorMessage(error));
      } finally {
        setBusy(false);
      }
    },
    [openFiles],
  );

  const onDropAccepted = useCallback((files: File[]) => void open(files), [open]);

  const {
    getRootProps,
    getInputProps,
    isDragActive,
    open: choose,
  } = useDropzone({
    noClick: true,
    noKeyboard: true,
    noDragEventsBubbling: true,
    onDropAccepted,
  });

  return (
    <div className="grid gap-3">
      <div
        {...getRootProps({
          role: "group",
          "aria-label": "Input files",
          className: cn(
            "border-input grid justify-items-center gap-3 rounded-lg border border-dashed px-6 py-10 text-center transition-colors",
            isDragActive && "border-primary bg-accent",
          ),
        })}
      >
        <input {...getInputProps()} aria-label="Input files to open" />
        <Upload aria-hidden className="text-muted-foreground size-6" />
        <p className="text-sm">
          Drop the tree, alignment and metadata files here, or a YAML configuration, or a run package.
        </p>
        <Button type="button" onClick={choose} disabled={busy}>
          {busy && <Spinner />}
          Choose files
        </Button>
        <SupportedFormats />
      </div>
      {failure !== undefined && (
        <Alert variant="destructive">
          <AlertDescription>The files cannot be opened: {failure}</AlertDescription>
        </Alert>
      )}
      {unplaced.length > 0 && (
        <Alert>
          <AlertTitle>
            {unplaced.length === 1 ? "One file was not added" : `${unplaced.length} files were not added`}
          </AlertTitle>
          <AlertDescription>
            <ul className="grid gap-0.5">
              {withListKeys(unplaced, ({ name }) => name).map(({ key, item: { name, reason } }) => (
                <li key={key}>
                  <span className="font-mono text-xs">{name}</span>: {reason}
                </li>
              ))}
            </ul>
          </AlertDescription>
        </Alert>
      )}
    </div>
  );
}

function SupportedFormats() {
  const slots = commandSettings(START_COMMAND).inputs;

  return (
    <dl className="text-muted-foreground grid grid-cols-[auto_auto] gap-x-3 gap-y-0.5 text-left text-xs">
      {slots.map((slot) => (
        <FormatRow key={slot.kind} label={slot.label} formats={slot.formats} />
      ))}
      <FormatRow label="Configuration" formats="YAML (.yaml, .yml)" />
      <FormatRow label="Run package" formats="ZIP (.zip), downloaded from a run" />
    </dl>
  );
}

function FormatRow({ label, formats }: { label: string; formats: string }) {
  return (
    <>
      <dt className="text-foreground">{label}</dt>
      <dd>{formats}</dd>
    </>
  );
}

function RecentWork() {
  const { all } = useListedRuns();
  const latest = useMemo(() => groupRuns(all).slice(0, RECENT_GROUPS), [all]);

  if (latest.length === 0) {
    return null;
  }

  return (
    <section aria-labelledby="recent-title" className="grid gap-3">
      <h2 id="recent-title" className="text-lg font-bold">
        Recent work
      </h2>
      <ul className="grid gap-1.5">
        {latest.map((group) => (
          <RecentGroup key={group.id} group={group} />
        ))}
      </ul>
    </section>
  );
}

function RecentGroup({ group }: { group: RunGroup }) {
  const newest = newestRun(group);

  if (newest === undefined) {
    return null;
  }

  const variants = group.entries.length - 1;

  return (
    <li>
      <Link
        to="/runs/$id/results"
        params={{ id: newest.id }}
        className="hover:bg-accent flex items-center gap-3 rounded-md border px-3 py-2 text-sm"
      >
        <span title={statusLabel(newest.status)}>
          <StatusIcon status={newest.status} />
        </span>
        <span className="min-w-0 flex-1 truncate font-bold">{newest.title}</span>
        {variants > 0 && (
          <span className="text-muted-foreground text-xs">
            {variants === 1 ? "1 more run" : `${variants} more runs`} in this line of work
          </span>
        )}
        <span className="text-muted-foreground text-xs tabular-nums">
          {newest.status === "ok" ? headlineText(newest.headline) : statusLabel(newest.status)}
        </span>
      </Link>
    </li>
  );
}
