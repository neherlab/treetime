import { InputPanel } from "./InputPanel";
import { ResultsPanel } from "./ResultsPanel";
import { useAppStore } from "../store/app-store";

export function Workspace() {
  const showResults = useAppStore((s) => s.showResults);
  const runStatus = useAppStore((s) => s.runStatus);

  return (
    <div className="flex flex-1 flex-col overflow-hidden lg:flex-row">
      <div className="shrink-0 overflow-y-auto border-b border-gray-200 lg:w-96 lg:border-b-0 lg:border-r dark:border-gray-700">
        <InputPanel />
      </div>

      <div className="flex flex-1 flex-col overflow-hidden">
        {showResults || runStatus === "completed" ? (
          <ResultsPanel />
        ) : (
          <div className="flex flex-1 items-center justify-center text-sm text-gray-400 dark:text-gray-600">
            <div className="text-center">
              <p className="mb-1">Results will appear here</p>
              <p className="text-xs">Load input files and click Run</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
