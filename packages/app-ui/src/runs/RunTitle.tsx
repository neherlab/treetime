import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsUpdate } from "@neherlab/app-contracts/client";
import { Pencil } from "lucide-react";
import { useCallback, useEffect, useRef, useState } from "react";

import { useApiMutation } from "../api/hooks";
import { Button, Toast } from "../ui";

export function RunTitle({ record }: { record: RunRecord }) {
  const toasts = Toast.useToastManager();
  const [draft, setDraft] = useState<string | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  const settled = useRef(true);
  const editing = draft !== null;

  useEffect(() => {
    if (editing) {
      inputRef.current?.focus();
      inputRef.current?.select();
    }
  }, [editing]);

  const { mutateAsync: rename } = useApiMutation((context, title: string) =>
    runsUpdate({ ...context, path: { id: record.id }, body: { title } }),
  );

  const startEditing = useCallback(() => {
    settled.current = false;
    setDraft(record.title);
  }, [record.title]);

  const cancel = useCallback(() => {
    settled.current = true;
    setDraft(null);
  }, []);

  const save = useCallback(async () => {
    if (settled.current) {
      return;
    }

    settled.current = true;
    const title = draft?.trim() ?? "";
    setDraft(null);

    if (title === "" || title === record.title) {
      return;
    }

    try {
      await rename(title);
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be renamed", description: errorMessage(error) });
    }
  }, [draft, record.title, rename, toasts]);

  const onChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => setDraft(event.target.value), []);

  const onBlur = useCallback(() => void save(), [save]);

  const onKeyDown = useCallback(
    (event: React.KeyboardEvent<HTMLInputElement>) => {
      if (event.key === "Enter") {
        event.preventDefault();
        void save();
      } else if (event.key === "Escape") {
        event.preventDefault();
        cancel();
      }
    },
    [cancel, save],
  );

  if (editing) {
    return (
      <input
        type="text"
        ref={inputRef}
        aria-label="Run title"
        value={draft}
        onChange={onChange}
        onBlur={onBlur}
        onKeyDown={onKeyDown}
        className="border-line-strong bg-surface-1 w-full rounded-md border px-2 text-2xl leading-tight font-bold"
      />
    );
  }

  return (
    <div className="flex items-center gap-1.5">
      <h1 className="text-2xl leading-tight font-bold">{record.title}</h1>
      <Button type="button" variant="ghost" size="sm" onClick={startEditing} aria-label="Rename the run" title="Rename">
        <Pencil size={14} aria-hidden />
      </Button>
    </div>
  );
}
