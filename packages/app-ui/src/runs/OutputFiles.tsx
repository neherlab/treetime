import { type RunFile, errorMessage, type OutputSelection } from "@neherlab/app-contracts";
import { runsFiles } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useCallback, useMemo } from "react";
import Download from "~icons/lucide/download";

import { useApiContext } from "../api/context";
import { useApi } from "../api/hooks";
import { formatBytes } from "../format";
import { useHost } from "../host-context";
import { browserFilePicker, downloadUrl, saveRunOutput, setBrowserTimer, type RunOutput } from "./saveRunOutput";

export function OutputFilesTable({ runId, files }: { runId: string; files: readonly RunFile[] }) {
  const perAlignment = files.some((file) => (file.alignment?.number ?? 0) > 1);

  return (
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
          <FileRow key={file.path} runId={runId} file={file} perAlignment={perAlignment} />
        ))}
      </TableBody>
    </Table>
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
  const pickFile = useMemo(() => browserFilePicker(), []);

  return useCallback(
    async (output: RunOutput, failureTitle: string) => {
      try {
        const reply = await saveRunOutput(
          { host, client, pickFile, download: downloadUrl, objectUrls: URL, setTimer: setBrowserTimer },
          output,
        );

        if (reply.kind === "error") {
          toasts.add({ title: failureTitle, description: [reply.error.message, ...reply.error.causes].join(": ") });
        }
      } catch (failure: unknown) {
        toasts.add({ title: failureTitle, description: errorMessage(failure) });
      }
    },
    [client, host, pickFile, toasts],
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
