import { Empty, EmptyDescription, EmptyHeader, EmptyTitle } from "@neherlab/design/ui/empty";

import { commandSettings } from "../settings/catalog";
import { inputFactsRequest } from "../settings/inputs";
import { useWorkspace } from "./workspaceContext";

export function DraftMain() {
  const { command, settled } = useWorkspace();
  const hasData = inputFactsRequest(command, settled) !== null;
  const { title } = commandSettings(command);

  return (
    <Empty className="h-full">
      <EmptyHeader>
        <EmptyTitle>{hasData ? `${title} not run yet` : "No data yet"}</EmptyTitle>
        <EmptyDescription>
          {hasData
            ? "Check the options on the left, then run the analysis. The results appear here."
            : "Add the input files in the Data section of the options, or open an example dataset."}
        </EmptyDescription>
      </EmptyHeader>
    </Empty>
  );
}
