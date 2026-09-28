import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord } from "@neherlab/app-contracts";
import { runsFiles } from "@neherlab/app-contracts/client";
import { Download } from "lucide-react";
import { useCallback, useState } from "react";

import { useApiContext } from "../api/context";
import { useApi } from "../api/hooks";
import { Panel } from "../components/Panel";
import { formatBytes } from "../format";
import { downloadName, totalSize, type RunFileEntry } from "../results/files";
import type { Citation } from "../results/types";
import { Alert, AlertDescription } from "../ui/alert";
import { Button } from "../ui/button";
import { Spinner } from "../ui/spinner";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { useToastManager } from "../ui/toast";

export function OutputFiles({ record, citation }: { record: RunRecord; citation: Citation }) {
  const { save } = useApiContext();
  const toasts = useToastManager();

  const { data: files, error } = useApi((context) => runsFiles({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const [busy, setBusy] = useState(false);

  const downloadArchive = useCallback(async () => {
    setBusy(true);

    try {
      await save.saveRunArchive(record.id, downloadName(record.title, ".zip"));
    } catch (failure: unknown) {
      toasts.add({ title: "The archive cannot be downloaded", description: errorMessage(failure) });
    } finally {
      setBusy(false);
    }
  }, [record.id, record.title, save, toasts]);

  const onArchive = useCallback(() => void downloadArchive(), [downloadArchive]);

  return (
    <Panel
      title="Output files"
      hint={files === undefined ? undefined : `${files.length} files, ${formatBytes(totalSize(files))}`}
      actions={
        <Button type="button" variant="outline" size="sm" onClick={onArchive} disabled={busy || files === undefined}>
          {busy ? <Spinner /> : <Download aria-hidden />}
          Download all (.zip)
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
              <TableHead>Contents</TableHead>
              <TableHead className="text-right">Size</TableHead>
              <TableHead>
                <span className="sr-only">Download</span>
              </TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {files.map((file) => (
              <FileRow key={file.path} runId={record.id} file={file} />
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

function FileRow({ runId, file }: { runId: string; file: RunFileEntry }) {
  const { save } = useApiContext();
  const toasts = useToastManager();

  const download = useCallback(async () => {
    try {
      await save.saveRunFile(runId, file.path, file.path.split("/").at(-1) ?? file.path);
    } catch (failure: unknown) {
      toasts.add({ title: "The file cannot be downloaded", description: errorMessage(failure) });
    }
  }, [file.path, runId, save, toasts]);

  const onDownload = useCallback(() => void download(), [download]);

  return (
    <TableRow>
      <TableCell>
        <code className="font-mono text-xs">{file.path}</code>
      </TableCell>
      <TableCell className="text-muted-foreground text-xs whitespace-normal">{file.description}</TableCell>
      <TableCell className="text-right text-xs">{formatBytes(file.size)}</TableCell>
      <TableCell className="text-right">
        <Button type="button" variant="ghost" size="icon-sm" aria-label={`Download ${file.path}`} onClick={onDownload}>
          <Download aria-hidden />
        </Button>
      </TableCell>
    </TableRow>
  );
}
