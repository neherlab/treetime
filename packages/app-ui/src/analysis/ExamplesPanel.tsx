import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand, Dataset, ExampleConfig, RunRecord } from "@neherlab/app-contracts";
import { datasets, runsGet, runsList } from "@neherlab/app-contracts/client";
import { useCallback, useMemo, useState } from "react";
import Database from "~icons/lucide/database";
import FileCog from "~icons/lucide/file-cog";

import { useApi, useApiQueries } from "../api/hooks";
import type { ApiCallContext } from "../api/keys";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { datasetInputs, runInputAssignments, type InputAssignment } from "../settings/inputs";
import { Alert, AlertDescription } from "../ui/alert";
import { Button } from "../ui/button";
import { Card } from "../ui/card";
import { Spinner } from "../ui/spinner";
import { Table, TableBody, TableCell, TableRow } from "../ui/table";
import { useToastManager } from "../ui/toast";
import { useConfigLoader } from "./useConfigLoader";
import { useInputActions } from "./useInputActions";

const EARLIER_RUNS_SHOWN = 8;

export const EXAMPLES_PANELS = {
  datasets: { title: "Example datasets", hint: "Bundled datasets with tree, alignment and metadata", icon: Database },
  configs: { title: "Example configs", hint: "Ready-made settings and inputs from the repository", icon: FileCog },
} as const;

export type ExamplesPanelKind = keyof typeof EXAMPLES_PANELS;

export const EXAMPLES_PANEL_KINDS: readonly ExamplesPanelKind[] = ["datasets", "configs"];

export function ExamplesPanel({
  kind,
  command,
  close,
}: {
  kind: ExamplesPanelKind;
  command: AppCommand;
  close: () => void;
}) {
  const { data: catalog, error } = useApi((context) => datasets(context), { staleTime: Infinity });

  return (
    <Card size="sm" className="gap-0 py-0">
      <PanelHeading title={EXAMPLES_PANELS[kind].title} hint={EXAMPLES_PANELS[kind].hint}>
        <Button type="button" variant="ghost" size="sm" onClick={close}>
          Close
        </Button>
      </PanelHeading>
      {error !== null && (
        <Alert variant="destructive" className="m-3.5 w-auto">
          <AlertDescription>The examples cannot be listed: {error.message}</AlertDescription>
        </Alert>
      )}
      {kind === "datasets" ? (
        <>
          <div className="max-h-72 overflow-auto overscroll-contain">
            <Table>
              <TableBody>
                {(catalog?.datasets ?? []).map((dataset) => (
                  <DatasetRow key={dataset.name} command={command} dataset={dataset} close={close} />
                ))}
              </TableBody>
            </Table>
          </div>
          <EarlierInputs command={command} close={close} />
        </>
      ) : (
        <div className="max-h-72 overflow-auto overscroll-contain">
          <Table>
            <TableBody>
              {(catalog?.examples ?? []).map((example) => (
                <ExampleRow key={example.path} example={example} />
              ))}
            </TableBody>
          </Table>
        </div>
      )}
    </Card>
  );
}

function DatasetRow({ command, dataset, close }: { command: AppCommand; dataset: Dataset; close: () => void }) {
  const actions = useInputActions(command);
  const inputs = useMemo(() => datasetInputs(dataset, command), [command, dataset]);

  const load = useCallback(() => {
    actions.assignAll(inputs, "dataset");
    close();
  }, [actions, close, inputs]);

  return (
    <TableRow>
      <TableCell className="font-mono text-xs">{dataset.name}</TableCell>
      <TableCell className="text-muted-foreground text-xs whitespace-normal">
        {inputs.map((input) => input.key).join(", ")}
      </TableCell>
      <TableCell className="text-right">
        <Button type="button" variant="outline" size="sm" disabled={inputs.length === 0} onClick={load}>
          Load
        </Button>
      </TableCell>
    </TableRow>
  );
}

function ExampleRow({ example }: { example: ExampleConfig }) {
  const loadConfig = useConfigLoader();
  const toasts = useToastManager();
  const [loading, setLoading] = useState(false);

  const load = useCallback(async () => {
    setLoading(true);

    try {
      const result = await loadConfig(example.content, example.command, false, example.folder);

      if (!result.loaded) {
        toasts.add({ title: `${example.path} cannot be loaded`, description: result.messages.join("; ") });
      }
    } catch (error: unknown) {
      toasts.add({
        title: `${example.path} cannot be loaded`,
        description: errorMessage(error),
      });
    } finally {
      setLoading(false);
    }
  }, [example, loadConfig, toasts]);

  const onLoad = useCallback(() => void load(), [load]);

  return (
    <TableRow>
      <TableCell className="whitespace-normal">{example.title}</TableCell>
      <TableCell className="text-muted-foreground font-mono text-xs whitespace-normal">{example.command}</TableCell>
      <TableCell className="text-right">
        <Button type="button" variant="outline" size="sm" disabled={loading} onClick={onLoad}>
          {loading && <Spinner />}
          Load
        </Button>
      </TableCell>
    </TableRow>
  );
}

function EarlierInputs({ command, close }: { command: AppCommand; close: () => void }) {
  const { data: runList } = useApi((context) => runsList(context));

  const ids = (runList?.runs ?? [])
    .filter((run) => run.status === "ok")
    .slice(0, EARLIER_RUNS_SHOWN)
    .map((run) => run.id);

  const records = useApiQueries(ids.map((id) => (context: ApiCallContext) => runsGet({ ...context, path: { id } })));
  const slots = new Set<string>(COMMAND_SETTINGS[command].inputs.map((input) => input.kind));

  const rows = records.flatMap((query) => {
    const record = query.data;

    if (record === undefined) {
      return [];
    }

    const inputs = runInputAssignments(command, record.inputs).filter((input) => slots.has(input.key));

    return inputs.length > 0 ? [{ record, inputs }] : [];
  });

  if (rows.length === 0) {
    return null;
  }

  return (
    <>
      <PanelHeading title="Inputs of earlier runs" hint="The files an earlier run read" />
      <Table>
        <TableBody>
          {rows.map(({ record, inputs }) => (
            <EarlierRow key={record.id} command={command} record={record} inputs={inputs} close={close} />
          ))}
        </TableBody>
      </Table>
    </>
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
    actions.assignAll(inputs, "run");
    close();
  }, [actions, close, inputs]);

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
