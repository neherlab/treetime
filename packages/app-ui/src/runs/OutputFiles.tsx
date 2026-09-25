import type { RunRecordResult } from "@neherlab/app-contracts";
import { Download } from "lucide-react";
import { useCallback, useState } from "react";

import { useBridge } from "../BridgeContext";
import { formatBytes } from "../format";
import { useRunFiles } from "../queries";
import { downloadName, fileDescription, totalSize, type RunFileEntry } from "../results/files";
import { CITATION, CITATION_DOI } from "../results/methods";
import { Button, Toast } from "../ui";
import { saveBytes } from "./download";
import { Panel } from "./Panel";
import { useCopy } from "./useCopy";

export function OutputFiles({ record, methods }: { record: RunRecordResult; methods: string | undefined }) {
  const bridge = useBridge();
  const toasts = Toast.useToastManager();
  const copy = useCopy();
  const { data: files, error } = useRunFiles(record.id, true);
  const [busy, setBusy] = useState(false);

  const downloadArchive = useCallback(async () => {
    setBusy(true);

    try {
      saveBytes(await bridge.runArchive(record.id), downloadName(record.title, ".zip"), "application/zip");
    } catch (failure: unknown) {
      toasts.add({
        title: "The archive cannot be downloaded",
        description: failure instanceof Error ? failure.message : "",
      });
    } finally {
      setBusy(false);
    }
  }, [bridge, record.id, record.title, toasts]);

  const onArchive = useCallback(() => void downloadArchive(), [downloadArchive]);

  const copyMethods = useCallback(() => {
    if (methods !== undefined) {
      copy(methods, "Methods text copied to the clipboard");
    }
  }, [copy, methods]);

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
      <div className="border-line text-ink-muted grid gap-2 border-t px-3.5 py-3 text-xs">
        <p className="m-0">
          Please cite: {CITATION}{" "}
          <a href={CITATION_DOI} target="_blank" rel="noopener noreferrer" className="text-accent font-bold">
            doi:10.1093/ve/vex042
          </a>
        </p>
        {methods !== undefined && (
          <div className="grid gap-1.5">
            <p className="text-ink m-0 leading-relaxed">{methods}</p>
            <div>
              <Button type="button" variant="ghost" size="sm" onClick={copyMethods}>
                Copy methods text
              </Button>
            </div>
          </div>
        )}
      </div>
    </Panel>
  );
}

function FileRow({ runId, file }: { runId: string; file: RunFileEntry }) {
  const bridge = useBridge();
  const toasts = Toast.useToastManager();

  const download = useCallback(async () => {
    try {
      saveBytes(
        await bridge.readRunFile(runId, file.path),
        file.path.split("/").at(-1) ?? file.path,
        "application/octet-stream",
      );
    } catch (failure: unknown) {
      toasts.add({
        title: "The file cannot be downloaded",
        description: failure instanceof Error ? failure.message : "",
      });
    }
  }, [bridge, file.path, runId, toasts]);

  const onDownload = useCallback(() => void download(), [download]);

  return (
    <tr className="border-line border-t">
      <td className="px-3.5 py-1.5">
        <code className="font-mono text-xs">{file.path}</code>
      </td>
      <td className="text-ink-muted px-3.5 py-1.5 text-xs">{fileDescription(file)}</td>
      <td className="px-3.5 py-1.5 text-right text-xs tabular-nums">{formatBytes(file.size)}</td>
      <td className="px-2 py-1 text-right">
        <Button type="button" variant="ghost" size="icon" aria-label={`Download ${file.path}`} onClick={onDownload}>
          <Download size={13} aria-hidden />
        </Button>
      </td>
    </tr>
  );
}
