import type { AppCommand, RunRecord } from "@neherlab/app-contracts";
import { runsGet } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Card } from "@neherlab/design/ui/card";
import { Table, TableBody, TableCell, TableRow } from "@neherlab/design/ui/table";
import { useCallback } from "react";

import { useApiQueries } from "../api/hooks";
import type { ApiCallContext } from "../api/keys";
import { commandSettings } from "../settings/catalog";
import { runInputAssignments, type InputAssignment } from "../settings/inputs";
import { useListedRuns } from "../shell/useListedRuns";
import { useInputActions } from "./useInputActions";

const EARLIER_RUNS_SHOWN = 8;

export function EarlierInputsPanel({ command, close }: { command: AppCommand; close: () => void }) {
  return (
    <Card size="sm" className="gap-0 py-0">
      <PanelHeading title="Inputs of earlier runs" hint="The files an earlier run read">
        <Button type="button" variant="ghost" size="sm" onClick={close}>
          Close
        </Button>
      </PanelHeading>
      <EarlierInputs command={command} close={close} />
    </Card>
  );
}

function EarlierInputs({ command, close }: { command: AppCommand; close: () => void }) {
  const { all } = useListedRuns();

  const ids = all
    .filter((run) => run.status === "ok")
    .slice(0, EARLIER_RUNS_SHOWN)
    .map((run) => run.id);

  const records = useApiQueries(ids.map((id) => (context: ApiCallContext) => runsGet({ ...context, path: { id } })));
  const slots = new Set<string>(commandSettings(command).inputs.map((input) => input.kind));

  const rows = records.flatMap((query) => {
    const record = query.data;

    if (record === undefined) {
      return [];
    }

    const inputs = runInputAssignments(command, record.inputs).filter((input) => slots.has(input.key));

    return inputs.length > 0 ? [{ record, inputs }] : [];
  });

  if (rows.length === 0) {
    return <p className="text-muted-foreground px-3.5 py-3 text-sm">No finished run read an input of this analysis.</p>;
  }

  return (
    <Table>
      <TableBody>
        {rows.map(({ record, inputs }) => (
          <EarlierRow key={record.id} command={command} record={record} inputs={inputs} close={close} />
        ))}
      </TableBody>
    </Table>
  );
}

function EarlierRow({
  command,
  record,
  inputs,
  close,
}: {
  command: AppCommand;
  record: RunRecord;
  inputs: readonly InputAssignment[];
  close: () => void;
}) {
  const actions = useInputActions(command);

  const use = useCallback(() => {
    actions.assignAll(inputs, { origin: { kind: "run" }, run_id: record.id });
    close();
  }, [actions, close, inputs, record.id]);

  return (
    <TableRow>
      <TableCell className="whitespace-normal">{record.title}</TableCell>
      <TableCell className="text-muted-foreground font-mono text-xs whitespace-normal">
        {inputs.map((input) => input.label).join(", ")}
      </TableCell>
      <TableCell className="text-right">
        <Button type="button" variant="outline" size="sm" onClick={use}>
          Use
        </Button>
      </TableCell>
    </TableRow>
  );
}

function PanelHeading({ title, hint, children }: { title: string; hint: string; children?: React.ReactNode }) {
  return (
    <div className="bg-muted/40 flex items-center gap-2.5 border-y px-3.5 py-2 first:border-t-0">
      <h3 className="text-sm font-bold">{title}</h3>
      <span className="text-muted-foreground text-xs">{hint}</span>
      {children !== undefined && <div className="ml-auto">{children}</div>}
    </div>
  );
}
