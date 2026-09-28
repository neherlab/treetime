import { errorMessage } from "@neherlab/app-contracts";
import type { RunRecord } from "@neherlab/app-contracts";
import { runsFiles } from "@neherlab/app-contracts/client";
import { Download } from "lucide-react";
import { useCallback, useState } from "react";

import { useApiContext } from "../api/context";
import { useApi } from "../api/hooks";
import { formatBytes } from "../format";
import { downloadName, totalSize, type RunFileEntry } from "../results/files";
import type { Citation } from "../results/types";
import { Button, Toast } from "../ui";
import { Panel } from "./Panel";

export function OutputFiles({ record, citation }: { record: RunRecord; citation: Citation }) {
  const { save } = useApiContext();
  const toasts = Toast.useToastManager();

  const { data: files, error } = useApi((context) => runsFiles({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const [busy, setBusy] = useState(false);

  const downloadArchive = useCallback(async () => {
    setBusy(true);

    try {
      await save.saveRunArchive(record.id, downloadName(record.title, ".zip"));
    } catch (failure: unknown) {
      toasts.add({
        title: "The archive cannot be downloaded",
        description: errorMessage(failure),
      });
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
          <Download size={13} aria-hidden />
          Download all (.zip)
        </Button>
      }
    >
      {error !== null && (
        <p className="text-signal-danger px-3.5 py-3">The file list cannot be loaded: {error.message}</p>
      )}
      {files !== undefined && (
        <table className="w-full border-collapse text-left">
          <thead>
            <tr className="text-ink-faint text-xs">
              <th className="px-3.5 py-1.5 font-normal">File</th>
              <th className="px-3.5 py-1.5 font-normal">Contents</th>
              <th className="px-3.5 py-1.5 text-right font-normal">Size</th>
              <th className="px-3.5 py-1.5" aria-label="Download" />
            </tr>
          </thead>
          <tbody>
            {files.map((file) => (
              <FileRow key={file.path} runId={record.id} file={file} />
            ))}
          </tbody>
        </table>
      )}
      <p className="border-line text-ink-muted m-0 border-t px-3.5 py-3 text-xs">
        Please cite: {citation.text}{" "}
        <a href={citation.url} target="_blank" rel="noopener noreferrer" className="text-accent font-bold">
          doi:{citation.doi}
        </a>
      </p>
    </Panel>
  );
}

function FileRow({ runId, file }: { runId: string; file: RunFileEntry }) {
  const { save } = useApiContext();
  const toasts = Toast.useToastManager();

  const download = useCallback(async () => {
    try {
      await save.saveRunFile(runId, file.path, file.path.split("/").at(-1) ?? file.path);
    } catch (failure: unknown) {
      toasts.add({
        title: "The file cannot be downloaded",
        description: errorMessage(failure),
      });
    }
  }, [file.path, runId, save, toasts]);

  const onDownload = useCallback(() => void download(), [download]);

  return (
    <tr className="border-line border-t">
      <td className="px-3.5 py-1.5">
        <code className="font-mono text-xs">{file.path}</code>
      </td>
      <td className="text-ink-muted px-3.5 py-1.5 text-xs">{file.description}</td>
      <td className="px-3.5 py-1.5 text-right text-xs tabular-nums">{formatBytes(file.size)}</td>
      <td className="px-2 py-1 text-right">
        <Button type="button" variant="ghost" size="icon" aria-label={`Download ${file.path}`} onClick={onDownload}>
          <Download size={13} aria-hidden />
        </Button>
      </td>
    </tr>
  );
}
