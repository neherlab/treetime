import { Alert, AlertDescription } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { CopyButton } from "@neherlab/design/ui/copy-button";
import { InputGroup, InputGroupAddon, InputGroupInput } from "@neherlab/design/ui/input-group";
import { useCallback, useMemo, useState } from "react";
import Download from "~icons/lucide/download";
import Search from "~icons/lucide/search";

import { OptionToggle } from "../components/OptionToggle";
import { countWarnings, filterLog, logText, type LogFilter, type RunProgress } from "../results/progress";
import { saveText } from "./saveRunOutput";
import { StageLog } from "./StageLog";

export function LogTab({
  progress,
  failure,
  fileName,
}: {
  progress: RunProgress;
  failure: string | undefined;
  fileName: string | undefined;
}) {
  const [filter, setFilter] = useState<LogFilter>("all");
  const [query, setQuery] = useState("");
  const entries = useMemo(() => filterLog(progress.entries, filter, query), [filter, progress.entries, query]);
  const lines = progress.entries.filter((entry) => entry.kind === "log").length;
  const warnings = countWarnings(progress.entries);

  const options = useMemo(
    () => [
      { value: "all" as const, label: `All ${lines}` },
      { value: "warnings" as const, label: `Warnings ${warnings}` },
    ],
    [lines, warnings],
  );

  const text = logText(entries);

  const onDownload = useCallback(() => {
    if (fileName !== undefined) {
      saveText(text, fileName, "text/plain");
    }
  }, [fileName, text]);

  const onQuery = useCallback((event: React.ChangeEvent<HTMLInputElement>) => setQuery(event.target.value), []);

  return (
    <StageLog
      progress={progress}
      filter={filter}
      query={query}
      hint="One card per stage"
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
          <CopyButton text={text} label="Copy the log" />
          <Button type="button" variant="outline" size="sm" onClick={onDownload} disabled={fileName === undefined}>
            <Download aria-hidden />
            Log (.txt)
          </Button>
        </>
      }
    >
      {failure !== undefined && (
        <Alert variant="destructive">
          <AlertDescription>The log cannot be followed: {failure}</AlertDescription>
        </Alert>
      )}
    </StageLog>
  );
}
