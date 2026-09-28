import { Search } from "lucide-react";
import { useCallback, useMemo, useState } from "react";

import { CopyButton } from "../components/CopyButton";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { countWarnings, filterLog, logText, type LogFilter, type RunProgress } from "../results/progress";
import { Alert, AlertDescription } from "../ui/alert";
import { InputGroup, InputGroupAddon, InputGroupInput } from "../ui/input-group";
import { LogLines } from "./LogLines";

export function LogTab({ progress, failure }: { progress: RunProgress; failure: string | undefined }) {
  const [filter, setFilter] = useState<LogFilter>("all");
  const [query, setQuery] = useState("");
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

  return (
    <Panel
      title="Log"
      actions={
        <>
          <OptionToggle label="Log filter" value={filter} onChange={setFilter} options={options} />
          <InputGroup className="h-8 w-56">
            <InputGroupAddon>
              <Search aria-hidden />
            </InputGroupAddon>
            <InputGroupInput
              type="search"
              value={query}
              onChange={onQuery}
              placeholder="Search the log"
              aria-label="Search the log"
            />
          </InputGroup>
          <CopyButton text={logText(entries)} label="Copy the log" />
        </>
      }
    >
      <div className="grid gap-2 p-2.5">
        {failure !== undefined && (
          <Alert variant="destructive">
            <AlertDescription>The log cannot be followed: {failure}</AlertDescription>
          </Alert>
        )}
        <LogLines
          entries={entries}
          className="h-[calc(100svh-var(--header-height)-16rem)] min-h-80"
          empty={progress.entries.length === 0 ? "No log lines yet." : "No matching lines."}
        />
      </div>
    </Panel>
  );
}
