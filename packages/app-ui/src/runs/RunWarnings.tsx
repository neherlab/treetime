import type { RunWarning } from "@neherlab/app-contracts";
import TriangleAlert from "~icons/lucide/triangle-alert";

import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Button } from "../ui/button";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";

export function RunWarnings({ warnings }: { warnings: readonly RunWarning[] }) {
  if (warnings.length === 0) {
    return null;
  }

  return (
    <Alert variant="warning">
      <TriangleAlert aria-hidden />
      <AlertTitle>The inputs have problems that make the results less reliable</AlertTitle>
      <AlertDescription>
        <ul className="grid gap-2">
          {warnings.map((warning) => (
            <li key={warning.message} className="grid gap-1">
              <p>{warning.message}</p>
              {warning.names.length > 0 && <WarningNames names={warning.names} />}
            </li>
          ))}
        </ul>
      </AlertDescription>
    </Alert>
  );
}

function WarningNames({ names }: { names: readonly string[] }) {
  return (
    <Collapsible className="grid gap-1">
      <CollapsibleTrigger
        render={<Button type="button" variant="link" size="sm" className="h-auto justify-self-start p-0" />}
      >
        {showNamesLabel(names.length)}
      </CollapsibleTrigger>
      <CollapsibleContent>
        <ul className="bg-muted/50 max-h-48 overflow-auto rounded-md border px-3 py-2 font-mono text-xs">
          {names.map((name) => (
            <li key={name}>{name}</li>
          ))}
        </ul>
      </CollapsibleContent>
    </Collapsible>
  );
}

export function showNamesLabel(count: number): string {
  return count === 1 ? "Show the name" : `Show the ${count} names`;
}
