import {
  type RunFile,
  errorMessage,
  type RunRecord,
  type Citation,
  type OutputSelection,
} from "@neherlab/app-contracts";
import { runsFiles } from "@neherlab/app-contracts/client";
import { Alert, AlertDescription } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { Spinner } from "@neherlab/design/ui/spinner";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useCallback, useState } from "react";
import Download from "~icons/lucide/download";

import { useApiContext } from "../api/context";
import { useApi } from "../api/hooks";
import { Panel } from "../components/Panel";
import { formatBytes } from "../format";
import { useHost } from "../host-context";
import { downloadName, totalSize } from "../results/files";
import { downloadUrl, saveRunOutput, setBrowserTimer, type RunOutput } from "./saveRunOutput";

export function OutputFiles({ record, citation }: { record: RunRecord; citation: Citation }) {
  const save = useSaveRunOutput();

  const { data: files, error } = useApi((context) => runsFiles({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const [busy, setBusy] = useState(false);

  const downloadArchive = useCallback(async () => {
    setBusy(true);
    await save({ id: record.id, name: downloadName(record.title, ".zip") }, "The run package cannot be downloaded");
    setBusy(false);
  }, [record.id, record.title, save]);

  const onArchive = useCallback(() => void downloadArchive(), [downloadArchive]);

  const perAlignment = files?.some((file) => (file.alignment?.number ?? 0) > 1) ?? false;

  return (
    <Panel
      title="Output files"
      hint={files === undefined ? undefined : `${files.length} files, ${formatBytes(totalSize(files))}`}
      actions={
        <Button type="button" variant="outline" size="sm" onClick={onArchive} disabled={busy || files === undefined}>
          {busy ? <Spinner /> : <Download aria-hidden />}
          Run package (.zip)
        </Button>
      }
    >
      {error !== null && (
        <Alert variant="destructive" className="m-3.5 w-auto">
          <AlertDescription>The file list cannot be loaded: {error.message}</AlertDescription>
        </Alert>
      )}
      {files !== undefined && (
        <Table>
          <TableHeader>
            <TableRow>
              <TableHead>File</TableHead>
              {perAlignment && <TableHead>Alignment</TableHead>}
              <TableHead>Contents</TableHead>
              <TableHead className="text-right">Size</TableHead>
              <TableHead>
                <span className="sr-only">Download</span>
              </TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {files.map((file) => (
              <FileRow key={file.path} runId={record.id} file={file} perAlignment={perAlignment} />
            ))}
          </TableBody>
        </Table>
      )}
      <p className="text-muted-foreground border-t px-3.5 py-3 text-xs">
        Please cite: {citation.text}{" "}
        <a
          href={citation.url}
          target="_blank"
          rel="noopener noreferrer"
          className="text-primary font-bold underline-offset-4 hover:underline"
        >
          doi:{citation.doi}
        </a>
      </p>
    </Panel>
  );
}

function FileRow({ runId, file, perAlignment }: { runId: string; file: RunFile; perAlignment: boolean }) {
  const save = useSaveRunOutput();

  const download = useCallback(async () => {
    await save({ id: runId, path: file.path, name: file.download_name }, "The file cannot be downloaded");
  }, [file.download_name, file.path, runId, save]);

  const onDownload = useCallback(() => void download(), [download]);

  return (
    <TableRow>
      <TableCell className="whitespace-normal">
        <span className={file.main ? "font-bold" : undefined}>{file.label}</span>
        {file.label !== file.path && <code className="text-muted-foreground block font-mono text-xs">{file.path}</code>}
      </TableCell>
      {perAlignment && (
        <TableCell className="text-xs">
          {file.alignment !== undefined && (
            <>
              {file.alignment.number} <code className="text-muted-foreground font-mono">{file.alignment.path}</code>
            </>
          )}
        </TableCell>
      )}
      <TableCell className="text-muted-foreground text-xs whitespace-normal">{file.description}</TableCell>
      <TableCell className="text-right text-xs">{formatBytes(file.size)}</TableCell>
      <TableCell className="text-right">
        <Button type="button" variant="ghost" size="icon-sm" aria-label={`Download ${file.label}`} onClick={onDownload}>
          <Download aria-hidden />
        </Button>
      </TableCell>
    </TableRow>
  );
}

export function useSaveRunOutput() {
  const { client } = useApiContext();
  const host = useHost();
  const toasts = useToastManager();

  return useCallback(
    async (output: RunOutput, failureTitle: string) => {
      try {
        const reply = await saveRunOutput(
          { host, client, download: downloadUrl, objectUrls: URL, setTimer: setBrowserTimer },
          output,
        );

        if (reply.kind === "error") {
          toasts.add({ title: failureTitle, description: [reply.error.message, ...reply.error.causes].join(": ") });
        }
      } catch (failure: unknown) {
        toasts.add({ title: failureTitle, description: errorMessage(failure) });
      }
    },
    [client, host, toasts],
  );
}

export function OutputDownload({
  runId,
  kind,
  children,
}: {
  runId: string;
  kind: OutputSelection;
  children: React.ReactNode;
}) {
  const save = useSaveRunOutput();

  const { data: files } = useApi((context) => runsFiles({ ...context, path: { id: runId } }), {
    staleTime: Infinity,
  });

  const file = files?.find((candidate) => candidate.kind === kind);

  const onDownload = useCallback(() => {
    if (file !== undefined) {
      void save({ id: runId, path: file.path, name: file.download_name }, "The file cannot be downloaded");
    }
  }, [file, runId, save]);

  if (file === undefined) {
    return null;
  }

  return (
    <Button type="button" variant="outline" size="sm" onClick={onDownload} title={file.label}>
      <Download aria-hidden />
      {children}
    </Button>
  );
}
