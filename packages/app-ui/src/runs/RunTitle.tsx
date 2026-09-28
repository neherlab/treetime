import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsUpdate } from "@neherlab/app-contracts/client";
import { Pencil } from "lucide-react";
import { useCallback, useRef, useState } from "react";

import { useApiMutation } from "../api/hooks";
import { Button } from "../ui/button";
import { Input } from "../ui/input";
import { useToastManager } from "../ui/toast";

export function RunTitle({ record }: { record: RunRecord }) {
  const toasts = useToastManager();
  const [draft, setDraft] = useState<string | null>(null);
  const settled = useRef(true);

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

  const focusInput = useCallback((input: HTMLInputElement | null) => {
    input?.focus();
    input?.select();
  }, []);

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

  if (draft !== null) {
    return (
      <Input
        type="text"
        aria-label="Run title"
        value={draft}
        onChange={onChange}
        ref={focusInput}
        onBlur={onBlur}
        onKeyDown={onKeyDown}
        className="font-heading h-9 text-xl font-semibold md:text-xl"
      />
    );
  }

  return (
    <div className="flex items-center gap-1">
      <h1 className="font-heading text-2xl leading-tight font-semibold">{record.title}</h1>
      <Button
        type="button"
        variant="ghost"
        size="icon-sm"
        onClick={startEditing}
        aria-label="Rename the run"
        title="Rename"
      >
        <Pencil aria-hidden />
      </Button>
    </div>
  );
}
