import { ArrowDown } from "lucide-react";
import { useCallback } from "react";
import { useStickToBottom } from "use-stick-to-bottom";

import { formatSeconds, type LogEntry } from "../results/progress";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";

const FOLLOW = { initial: "instant", resize: "instant" } as const;

const STILL = { initial: false, resize: "instant" } as const;

export function LogLines({
  entries,
  className,
  empty,
  follow,
}: {
  entries: readonly LogEntry[];
  className?: string;
  empty: string;
  follow: boolean;
}) {
  const { scrollRef, contentRef, isAtBottom, scrollToBottom } = useStickToBottom(follow ? FOLLOW : STILL);
  const jumpToEnd = useCallback(() => void scrollToBottom(), [scrollToBottom]);

  return (
    <div className={cn("relative flex min-h-0 flex-col", className)}>
      <div
        ref={scrollRef}
        className="bg-muted/50 min-h-0 flex-auto overflow-auto overscroll-contain rounded-md border font-mono text-xs leading-relaxed"
      >
        <div ref={contentRef} role="log" aria-live="off" className="px-3 py-2 wrap-anywhere whitespace-pre-wrap">
          {entries.length === 0 && <p className="text-muted-foreground">{empty}</p>}
          {entries.map((entry) => (
            <div key={`${entry.kind}:${entry.seq}`} className={lineClass(entry)}>
              <span className="text-muted-foreground select-none">{formatSeconds(entry.seconds)} </span>
              {entry.level === "warn" && <span className="font-bold">warning: </span>}
              {entry.level === "error" && <span>error: </span>}
              {entry.message}
            </div>
          ))}
        </div>
      </div>
      {follow && !isAtBottom && (
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
    return "text-primary font-bold";
  }

  if (entry.level === "warn") {
    return "text-warning";
  }

  if (entry.level === "error") {
    return "text-destructive font-bold";
  }

  return entry.level === "debug" || entry.level === "trace" ? "text-muted-foreground" : undefined;
}
