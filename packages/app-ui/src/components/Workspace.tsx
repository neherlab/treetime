import { useAppStore } from "../store/app-store";
import { InputPanel } from "./InputPanel";
import { ResultsPanel } from "./ResultsPanel";

export function Workspace() {
  const showResults = useAppStore((s) => s.showResults);
  const runStatus = useAppStore((s) => s.runStatus);

  return (
    <div className="flex flex-1 flex-col overflow-hidden lg:flex-row">
      <div className="border-line shrink-0 overflow-y-auto border-b lg:w-96 lg:border-r lg:border-b-0">
        <InputPanel />
      </div>

      <div className="flex flex-1 flex-col overflow-hidden">
        {showResults || runStatus === "completed" ? (
          <ResultsPanel />
        ) : (
          <div className="text-ink-faint flex flex-1 items-center justify-center text-sm">
            <div className="text-center">
              <p className="mb-1">Results will appear here</p>
              <p className="text-2xs">Load input files and start a run</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
