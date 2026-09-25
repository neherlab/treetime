import { Copy } from "lucide-react";
import { useCallback, useMemo, useState } from "react";

import { countWarnings, filterLog, logText, type LogFilter, type RunProgress } from "../results/progress";
import { Button, Segmented } from "../ui";
import { LogLines } from "./LogLines";
import { Panel } from "./Panel";
import { useCopy } from "./useCopy";

export function LogTab({ progress, failure }: { progress: RunProgress; failure: string | undefined }) {
  const [filter, setFilter] = useState<LogFilter>("all");
  const [query, setQuery] = useState("");
  const copy = useCopy();
  const entries = useMemo(() => filterLog(progress.entries, filter, query), [filter, progress.entries, query]);
  const warnings = countWarnings(progress.entries);
  const stages = progress.entries.filter((entry) => entry.kind === "stage").length;

  const options = useMemo(
    () => [
      { value: "all" as const, label: `All ${progress.entries.length}` },
      { value: "warnings" as const, label: `Warnings ${warnings}` },
      { value: "stages" as const, label: `Stages ${stages}` },
    ],
    [progress.entries.length, stages, warnings],
  );

  const onQuery = useCallback((event: React.ChangeEvent<HTMLInputElement>) => setQuery(event.target.value), []);
  const onCopy = useCallback(() => copy(logText(entries), "Log copied to the clipboard"), [copy, entries]);

  return (
    <Panel
      title="Log"
      actions={
        <>
          <Segmented label="Log filter" value={filter} onChange={setFilter} options={options} />
          <input
            type="search"
            value={query}
            onChange={onQuery}
            placeholder="Search the log"
            aria-label="Search the log"
            className="border-line-strong bg-surface-2 w-56 rounded-md border px-2 py-1 text-sm"
          />
          <Button type="button" variant="outline" size="sm" onClick={onCopy}>
            <Copy size={13} aria-hidden />
            Copy
          </Button>
        </>
      }
    >
      <div className="p-2.5">
        {failure !== undefined && <p className="text-signal-danger mb-2">The log cannot be followed: {failure}</p>}
        <LogLines
          entries={entries}
          className="max-h-[70vh]"
          empty={progress.entries.length === 0 ? "No log lines yet." : "No matching lines."}
        />
      </div>
    </Panel>
  );
}
