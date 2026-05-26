import { Upload, X } from "lucide-react";
import { clsx } from "clsx";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store/app-store";
import type { FileSlotConfig } from "../types";

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

  return (
    <div
      className={clsx(
        "rounded-lg border-2 border-dashed p-3 transition-colors",
        !relevant && "opacity-40",
        file
          ? "border-green-300 bg-green-50 dark:border-green-800 dark:bg-green-950/30"
          : "border-gray-200 bg-white hover:border-gray-300 dark:border-gray-700 dark:bg-gray-900 dark:hover:border-gray-600",
      )}
      onDrop={handleDrop}
      onDragOver={handleDragOver}
    >
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <span className="text-sm font-medium text-gray-700 dark:text-gray-300">{config.label}</span>
          {required && relevant && <span className="text-xs text-red-500">required</span>}
        </div>
        {file && (
          <button type="button" onClick={handleClear} className="rounded p-0.5 text-gray-400 hover:text-red-500">
            <X size={14} />
          </button>
        )}
      </div>

      {file ? (
        <div className="mt-1 flex items-center gap-2">
          <span className="truncate text-xs text-gray-600 dark:text-gray-400">{file.name}</span>
          <span className="shrink-0 text-xs text-gray-400">{formatSize(file.size)}</span>
        </div>
      ) : (
        <button
          type="button"
          onClick={() => inputRef.current?.click()}
          className="mt-1.5 flex w-full items-center justify-center gap-1.5 rounded-md bg-gray-50 py-2 text-xs text-gray-500 hover:bg-gray-100 dark:bg-gray-800 dark:text-gray-500 dark:hover:bg-gray-700"
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
