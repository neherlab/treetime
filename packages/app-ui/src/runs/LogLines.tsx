import { ArrowDown } from "lucide-react";
import { useCallback } from "react";
import { useStickToBottom } from "use-stick-to-bottom";

import { formatSeconds, type LogEntry } from "../results/progress";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";

const FOLLOW = { initial: "instant", resize: "instant" } as const;

export function LogLines({
  entries,
  className,
  empty,
}: {
  entries: readonly LogEntry[];
  className?: string;
  empty: string;
}) {
  const { scrollRef, contentRef, isAtBottom, scrollToBottom } = useStickToBottom(FOLLOW);
  const jumpToEnd = useCallback(() => void scrollToBottom(), [scrollToBottom]);

  return (
    <div className={cn("relative min-h-0", className)}>
      <div
        ref={scrollRef}
        className="bg-muted/50 size-full overflow-auto overscroll-contain rounded-md border font-mono text-xs leading-relaxed"
      >
        <div ref={contentRef} role="log" aria-live="off" className="px-3 py-2 wrap-anywhere whitespace-pre-wrap">
          {entries.length === 0 && <p className="text-muted-foreground">{empty}</p>}
          {entries.map((entry) => (
            <div key={`${entry.kind}:${entry.seq}`} className={lineClass(entry)}>
              <span className="text-muted-foreground select-none">{formatSeconds(entry.seconds)} </span>
              {entry.level === "warn" && <span className="font-medium">warning: </span>}
              {entry.level === "error" && <span>error: </span>}
              {entry.message}
            </div>
          ))}
        </div>
      </div>
      {!isAtBottom && (
        <Button type="button" size="xs" variant="secondary" onClick={jumpToEnd} className="absolute right-3 bottom-3">
          <ArrowDown aria-hidden />
          Follow new lines
        </Button>
      )}
    </div>
  );
}

function lineClass(entry: LogEntry): string | undefined {
  if (entry.kind === "stage") {
    return "text-primary font-medium";
  }

  if (entry.level === "warn") {
    return "text-warning";
  }

  if (entry.level === "error") {
    return "text-destructive font-medium";
  }

  return entry.level === "debug" || entry.level === "trace" ? "text-muted-foreground" : undefined;
}
