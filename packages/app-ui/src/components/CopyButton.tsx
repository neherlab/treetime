import { useClipboard } from "@mantine/hooks";
import { useCallback } from "react";
import Check from "~icons/lucide/check";
import Copy from "~icons/lucide/copy";

import { Button } from "../ui/button";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

export function CopyButton({ text, label, disabled }: { text: string; label: string; disabled?: boolean }) {
  const { copy, copied, error } = useClipboard();
  const onCopy = useCallback(() => copy(text), [copy, text]);
  const status = error === null ? (copied ? "Copied" : label) : "The browser blocked clipboard access";

  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button type="button" variant="ghost" size="icon-sm" aria-label={label} onClick={onCopy} disabled={disabled}>
            {copied ? <Check aria-hidden className="text-success" /> : <Copy aria-hidden />}
          </Button>
        }
      />
      <TooltipContent>{status}</TooltipContent>
    </Tooltip>
  );
}
