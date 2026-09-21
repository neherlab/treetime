import { Upload, X } from "lucide-react";
import { useCallback, useRef } from "react";

import { useAppStore } from "../store/app-store";
import type { FileSlotConfig } from "../types";
import { cn } from "../ui";

interface FileSlotProps {
  config: FileSlotConfig;
  relevant: boolean;
  required: boolean;
}

export function FileSlot({ config, relevant, required }: FileSlotProps) {
  const file = useAppStore((s) => s.files[config.kind]);
  const setFile = useAppStore((s) => s.setFile);
  const inputRef = useRef<HTMLInputElement>(null);

  const handleFileSelect = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const selected = e.target.files?.[0];

      if (selected) {
        setFile(config.kind, { name: selected.name, size: selected.size });
      }
    },
    [config.kind, setFile],
  );

  const handleClear = useCallback(() => {
    setFile(config.kind, undefined);

    if (inputRef.current) {
      inputRef.current.value = "";
    }
  }, [config.kind, setFile]);

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      const dropped = e.dataTransfer.files[0];

      if (dropped) {
        setFile(config.kind, { name: dropped.name, size: dropped.size });
      }
    },
    [config.kind, setFile],
  );

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
  }, []);

  const handleBrowse = useCallback(() => {
    inputRef.current?.click();
  }, []);

  return (
    <div
      className={cn(
        "rounded-md border border-dashed p-3 transition-colors",
        !relevant && "opacity-40",
        file ? "border-signal-ok/60 bg-signal-ok/10" : "border-line bg-surface-0 hover:border-line-strong",
      )}
      onDrop={handleDrop}
      onDragOver={handleDragOver}
    >
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <span className="text-ink text-sm font-medium">{config.label}</span>
          {required && relevant && <span className="text-2xs text-signal-danger">required</span>}
        </div>
        {file && (
          <button
            type="button"
            onClick={handleClear}
            aria-label={`Clear ${config.label}`}
            className="text-ink-faint hover:text-signal-danger rounded-sm p-0.5"
          >
            <X size={14} />
          </button>
        )}
      </div>

      {file ? (
        <div className="mt-1 flex items-center gap-2">
          <span className="text-2xs text-ink-muted truncate font-mono">{file.name}</span>
          <span className="text-2xs text-ink-faint shrink-0 font-mono">{formatSize(file.size)}</span>
        </div>
      ) : (
        <button
          type="button"
          onClick={handleBrowse}
          className="bg-surface-2 text-2xs text-ink-muted hover:bg-surface-3 mt-1.5 flex w-full items-center justify-center gap-1.5 rounded-md py-2"
        >
          <Upload size={14} />
          <span>{config.description}</span>
        </button>
      )}

      <input ref={inputRef} type="file" accept={config.accept} onChange={handleFileSelect} className="hidden" />
    </div>
  );
}

function formatSize(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;

  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;

  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}
