import { useCallback, useEffect, useRef } from "react";

import { formatSeconds, type LogEntry } from "../results/progress";
import { cn } from "../ui";

const STICK_THRESHOLD_PX = 8;

export function LogLines({
  entries,
  className,
  empty,
}: {
  entries: readonly LogEntry[];
  className?: string;
  empty: string;
}) {
  const box = useRef<HTMLDivElement>(null);
  const following = useRef(true);

  const onScroll = useCallback(() => {
    const element = box.current;

    if (element !== null) {
      following.current = element.scrollTop + element.clientHeight >= element.scrollHeight - STICK_THRESHOLD_PX;
    }
  }, []);

  useEffect(() => {
    const element = box.current;

    if (element !== null && following.current && entries.length > 0) {
      element.scrollTop = element.scrollHeight;
    }
  }, [entries]);

  return (
    <div
      ref={box}
      onScroll={onScroll}
      role="log"
      aria-live="off"
      className={cn(
        "border-line bg-surface-2 overflow-auto rounded-md border px-2.5 py-2 font-mono text-[0.71875rem] leading-relaxed whitespace-pre-wrap",
        className,
      )}
    >
      {entries.length === 0 && <p className="text-ink-faint m-0">{empty}</p>}
      {entries.map((entry) => (
        <div
          key={`${entry.kind}:${entry.seq}`}
          className={cn(
            entry.kind === "stage" && "text-accent font-bold",
            entry.level === "warn" && "text-signal-warn",
            entry.level === "error" && "text-signal-danger font-bold",
            (entry.level === "debug" || entry.level === "trace") && "text-ink-faint",
          )}
        >
          <span className="text-ink-faint select-none">{formatSeconds(entry.seconds)} </span>
          {entry.level === "warn" && <span className="font-bold">warning: </span>}
          {entry.level === "error" && <span>error: </span>}
          {entry.message}
        </div>
      ))}
    </div>
  );
}
