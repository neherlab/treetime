import type { LogLevel } from "@neherlab/app-contracts";
import { useVirtualizer } from "@tanstack/react-virtual";
import { useCallback } from "react";
import { useStickToBottom } from "use-stick-to-bottom";
import ArrowDown from "~icons/lucide/arrow-down";

import { formatSeconds, type LogEntry } from "../results/progress";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";

const LEVEL_CLASS: Record<LogLevel, string | undefined> = {
  trace: "text-muted-foreground",
  debug: "text-muted-foreground",
  info: undefined,
  warn: "text-warning",
  error: "text-destructive font-bold",
};

const FOLLOW = { initial: "instant", resize: "instant" } as const;

const STILL = { initial: false, resize: "instant" } as const;

const LINE_HEIGHT_ESTIMATE = 20;

const OVERSCAN = 20;

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
  const getScrollElement = useCallback(() => scrollRef.current, [scrollRef]);
  const estimateSize = useCallback(() => LINE_HEIGHT_ESTIMATE, []);
  const getItemKey = useCallback((index: number) => entryKey(entries[index]) ?? index, [entries]);

  // oxlint-disable-next-line react/incompatible-library -- no build enables React Compiler, and the virtualizer reaches only DOM refs, never a memoized child
  const virtualizer = useVirtualizer({
    count: entries.length,
    getScrollElement,
    estimateSize,
    getItemKey,
    overscan: OVERSCAN,
  });

  const rows = virtualizer.getVirtualItems();

  return (
    <div className={cn("relative flex min-h-0 flex-col", className)}>
      <div
        ref={scrollRef}
        className="bg-muted/50 min-h-0 flex-auto overflow-auto overscroll-contain rounded-md border font-mono text-xs leading-relaxed"
      >
        <div ref={contentRef} role="log" aria-live="off" className="px-3 py-2 wrap-anywhere whitespace-pre-wrap">
          {entries.length === 0 && <p className="text-muted-foreground">{empty}</p>}
          <div
            className="relative h-(--log-height)"
            // oxlint-disable-next-line react/forbid-dom-props -- the virtual list height and row offset are runtime values that no Tailwind class can hold
            style={{
              "--log-height": `${virtualizer.getTotalSize()}px`,
              "--log-offset": `${rows[0]?.start ?? 0}px`,
            }}
          >
            <div className="absolute inset-x-0 top-0 translate-y-(--log-offset)">
              {rows.map((row) => {
                const entry = entries[row.index];

                return (
                  entry !== undefined && (
                    <div
                      key={row.key}
                      data-index={row.index}
                      ref={virtualizer.measureElement}
                      className={lineClass(entry)}
                    >
                      <span className="text-muted-foreground select-none">{formatSeconds(entry.seconds)} </span>
                      {entry.level === "warn" && <span className="font-bold">warning: </span>}
                      {entry.level === "error" && <span>error: </span>}
                      {entry.message}
                    </div>
                  )
                );
              })}
            </div>
          </div>
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

function entryKey(entry: LogEntry | undefined): string | undefined {
  return entry === undefined ? undefined : `${entry.kind}:${entry.seq}`;
}

function lineClass(entry: LogEntry): string | undefined {
  return entry.kind === "stage" ? "text-primary font-bold" : LEVEL_CLASS[entry.level];
}

declare module "react" {
  interface CSSProperties {
    "--log-height"?: string;
    "--log-offset"?: string;
  }
}
