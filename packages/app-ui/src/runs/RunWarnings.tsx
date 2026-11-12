import type { RunWarning, RunWarningKind } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import TriangleAlert from "~icons/lucide/triangle-alert";

export type WarnedEstimate = "samples";

const WARNED_ESTIMATES: Readonly<Record<RunWarningKind, WarnedEstimate | undefined>> = {
  "duplicate-node-names": undefined,
  "duplicate-sequence-names": undefined,
  "duplicate-metadata-names": undefined,
  "duplicate-partition-names": undefined,
  "branch-lengths-not-converged": undefined,
  "dates-outside-input-dates": "samples",
  "filtered-tips-removed": "samples",
  "undated-tips-removed": "samples",
  "locations-without-coordinates": undefined,
};

export function RunWarnings({ warnings }: { warnings: readonly RunWarning[] }) {
  const summary = warningSummary(warnings);

  if (summary === undefined) {
    return null;
  }

  return (
    <Collapsible className="grid gap-2">
      <div className="text-warning flex flex-wrap items-center gap-x-2 gap-y-1 text-sm">
        <TriangleAlert aria-hidden className="size-4 shrink-0" />
        <span className="min-w-0">{summary}</span>
        <CollapsibleTrigger render={<Button type="button" variant="link" size="sm" className="h-auto p-0" />}>
          Show
        </CollapsibleTrigger>
      </div>
      <CollapsibleContent>
        <ul className="bg-warning/5 border-warning/30 grid gap-2 rounded-md border px-3 py-2 text-sm">
          {warnings.map((warning) => (
            <li key={warning.message} className="grid gap-1">
              <p>{warning.message}</p>
              {warning.names.length > 0 && <WarningNames names={warning.names} />}
            </li>
          ))}
        </ul>
      </CollapsibleContent>
    </Collapsible>
  );
}

export function WarningNames({ names }: { names: readonly string[] }) {
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

export function warningSummary(warnings: readonly RunWarning[]): string | undefined {
  const first = warnings[0];

  if (first === undefined) {
    return undefined;
  }

  const count = warnings.length === 1 ? "1 warning" : `${warnings.length} warnings`;

  return `${count}: ${first.message.trim().replace(/\.$/u, "")}`;
}

export function warningsOf(warnings: readonly RunWarning[], estimate: WarnedEstimate): RunWarning[] {
  return warnings.filter((warning) => WARNED_ESTIMATES[warning.kind] === estimate);
}

export function showNamesLabel(count: number): string {
  return count === 1 ? "Show the name" : `Show the ${count} names`;
}
