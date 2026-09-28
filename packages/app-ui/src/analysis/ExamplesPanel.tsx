import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand, Dataset, ExampleConfig, RunRecord } from "@neherlab/app-contracts";
import { datasets, runsGet, runsList } from "@neherlab/app-contracts/client";
import { useCallback, useMemo, useState } from "react";

import { useApi, useApiQueries } from "../api/hooks";
import type { ApiCallContext } from "../api/keys";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { datasetInputs, runInputAssignments, type InputAssignment } from "../settings/inputs";
import { Button, Toast } from "../ui";
import { useConfigLoader } from "./useConfigLoader";
import { useInputActions } from "./useInputActions";

const EARLIER_RUNS_SHOWN = 8;

export function ExamplesPanel({ command, close }: { command: AppCommand; close: () => void }) {
  const { data: catalog, error } = useApi((context) => datasets(context), { staleTime: Infinity });

  return (
    <div className="border-line bg-surface-1 rounded-lg border">
      <PanelHeading title="Example data" hint="Bundled datasets with tree, alignment and metadata">
        <Button type="button" variant="ghost" size="sm" onClick={close}>
          Close
        </Button>
      </PanelHeading>
      {error !== null && (
        <p className="text-signal-danger px-3.5 py-2">The examples cannot be listed: {error.message}</p>
      )}
      <div className="max-h-72 overflow-auto">
        <table className="w-full border-collapse">
          <tbody>
            {(catalog?.datasets ?? []).map((dataset) => (
              <DatasetRow key={dataset.name} command={command} dataset={dataset} close={close} />
            ))}
          </tbody>
        </table>
      </div>
      <PanelHeading title="Example configs" hint="Ready-made settings and inputs from the repository" />
      <table className="w-full border-collapse">
        <tbody>
          {(catalog?.examples ?? []).map((example) => (
            <ExampleRow key={example.path} example={example} />
          ))}
        </tbody>
      </table>
      <EarlierInputs command={command} close={close} />
    </div>
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
    <tr className="border-line border-b">
      <td className="px-3.5 py-1.5 font-mono text-xs">{dataset.name}</td>
      <td className="text-ink-faint px-2 py-1.5 text-xs">{inputs.map((input) => input.key).join(", ")}</td>
      <td className="px-3.5 py-1 text-right">
        <Button type="button" variant="outline" size="sm" disabled={inputs.length === 0} onClick={load}>
          Load
        </Button>
      </td>
    </tr>
  );
}

function ExampleRow({ example }: { example: ExampleConfig }) {
  const loadConfig = useConfigLoader();
  const toasts = Toast.useToastManager();
  const [loading, setLoading] = useState(false);

  const load = useCallback(async () => {
    setLoading(true);

    try {
      const result = await loadConfig(example.content, example.command, false);

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
    <tr className="border-line border-b">
      <td className="px-3.5 py-1.5">{example.title}</td>
      <td className="text-ink-faint px-2 py-1.5 font-mono text-xs">{example.command}</td>
      <td className="px-3.5 py-1 text-right">
        <Button type="button" variant="outline" size="sm" disabled={loading} onClick={onLoad}>
          {loading ? "Loading..." : "Load"}
        </Button>
      </td>
    </tr>
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
      <table className="w-full border-collapse">
        <tbody>
          {rows.map(({ record, inputs }) => (
            <EarlierRow key={record.id} command={command} record={record} inputs={inputs} close={close} />
          ))}
        </tbody>
      </table>
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
    <tr className="border-line border-b">
      <td className="px-3.5 py-1.5">{record.title}</td>
      <td className="text-ink-faint px-2 py-1.5 font-mono text-xs">{inputs.map((input) => input.label).join(", ")}</td>
      <td className="px-3.5 py-1 text-right">
        <Button type="button" variant="outline" size="sm" onClick={use}>
          Use
        </Button>
      </td>
    </tr>
  );
}

function PanelHeading({ title, hint, children }: { title: string; hint: string; children?: React.ReactNode }) {
  return (
    <div className="border-line flex items-center gap-2.5 border-b px-3.5 py-2.5">
      <h3 className="font-bold">{title}</h3>
      <span className="text-ink-faint text-xs">{hint}</span>
      {children !== undefined && <div className="ml-auto">{children}</div>}
    </div>
  );
}
