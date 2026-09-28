import { useLayoutEffect, useRef } from "react";

import { formatSeconds, type LogEntry } from "../results/progress";
import { cn } from "../ui";
import { isScrolledToEnd } from "./logFollow";

export type LogScroller = "box" | "page";

const SCROLLERS: Record<LogScroller, Scroller> = {
  box: { className: "overflow-auto", element: (box) => box, followsFirstLines: () => true },
  page: { className: undefined, element: () => document.scrollingElement, followsFirstLines: isScrolledToEnd },
};

export function LogLines({
  entries,
  scroller,
  className,
  empty,
}: {
  entries: readonly LogEntry[];
  scroller: LogScroller;
  className?: string;
  empty: string;
}) {
  const box = useRef<HTMLDivElement>(null);
  const renderedHeight = useRef<number | undefined>(undefined);
  const { className: scrollerClassName, element, followsFirstLines } = SCROLLERS[scroller];

  useLayoutEffect(() => {
    const node = box.current;
    const scrolled = node === null ? null : element(node);

    if (scrolled === null || entries.length === 0) {
      return;
    }

    const previousHeight = renderedHeight.current;

    const follows =
      previousHeight === undefined
        ? followsFirstLines(scrolled)
        : isScrolledToEnd({
            scrollTop: scrolled.scrollTop,
            clientHeight: scrolled.clientHeight,
            scrollHeight: previousHeight,
          });

    if (follows) {
      scrolled.scrollTop = scrolled.scrollHeight;
    }

    renderedHeight.current = scrolled.scrollHeight;
  }, [element, entries, followsFirstLines]);

  return (
    <div
      ref={box}
      role="log"
      aria-live="off"
      className={cn(
        "border-line bg-surface-2 rounded-md border px-2.5 py-2 font-mono text-[0.71875rem] leading-relaxed wrap-anywhere whitespace-pre-wrap",
        scrollerClassName,
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

interface Scroller {
  className: string | undefined;
  element: (box: HTMLDivElement) => Element | null;
  followsFirstLines: (scrolled: Element) => boolean;
}
