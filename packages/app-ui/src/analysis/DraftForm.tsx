import type { AppCommand, JobId } from "@neherlab/app-contracts";
import { configCheck, inputsCheck, runConfig as resolveRunConfig, runsList } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { useHotkey } from "@tanstack/react-hotkeys";
import { useDebouncedValue } from "@tanstack/react-pacer";
import { keepPreviousData } from "@tanstack/react-query";
import { Outlet } from "@tanstack/react-router";
import { useCallback, useEffect, useMemo, useState } from "react";
import { FormProvider, useForm, useWatch } from "react-hook-form";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { useApi } from "../api/hooks";
import { RUN_HOTKEY } from "../hotkeys";
import { commandSettings } from "../settings/catalog";
import { COMMAND_VERBS } from "../settings/commands";
import { normalizeConfig } from "../settings/config";
import { inputFactsRequest } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { RunFooter } from "./ChecksPanel";
import { configResolver } from "./configResolver";
import { datasetLabel } from "./datasetLabel";
import { fromFormConfig, toFormConfig, type FormConfig } from "./formValues";
import { OptionsPanel } from "./OptionsPanel";
import { useStartRun } from "./useStartRun";
import { WorkspaceContext, type WorkspaceView } from "./workspaceContext";
import { WorkspaceLayout } from "./WorkspaceLayout";

const CHECK_DEBOUNCE = { wait: 400 };

export function DraftForm({ command, shownRunId }: { command: AppCommand; shownRunId: JobId | undefined }) {
  const specs = commandSettings(command).settings;
  const [initial] = useState(() => normalizeConfig(specs, useDraftStore.getState().draft.config));

  const form = useForm<FormConfig>({
    defaultValues: toFormConfig(initial),
    resolver: configResolver(command),
    mode: "onChange",
  });

  const setConfig = useDraftStore((state) => state.setConfig);
  const readOnly = useDraftStore((state) => state.viewing !== undefined);
  const watched = useWatch({ control: form.control });

  const config = useMemo(() => normalizeConfig(specs, fromFormConfig(watched)), [specs, watched]);

  const [settled] = useDebouncedValue(config, CHECK_DEBOUNCE);
  const factsRequest = useMemo(() => inputFactsRequest(command, settled), [command, settled]);

  const { data: facts } = useApi(
    (context) => inputsCheck({ ...context, body: factsRequest ?? { command, config: {} } }),
    { enabled: factsRequest !== null, placeholderData: keepPreviousData, staleTime: 30_000 },
  );

  const choices = useDraftStore((state) => state.draft.choices);

  const { data: check, isPlaceholderData: checkPending } = useApi(
    (context) => {
      const text = JSON.stringify(settled);
      const picked = choices ?? [];

      return configCheck({
        ...context,
        body:
          facts === undefined
            ? { command, text, choices: picked }
            : { command, text, choices: picked, input_facts: facts },
      });
    },
    { placeholderData: keepPreviousData, staleTime: Infinity },
  );

  const { data: runConfig } = useApi(
    (context) => resolveRunConfig({ ...context, body: { command, config: settled } }),
    { placeholderData: keepPreviousData, staleTime: 10_000 },
  );

  const { data: runList } = useApi((context) => runsList(context));
  const startRun = useStartRun(command);

  useEffect(
    () =>
      form.subscribe({
        formState: { values: true },
        callback: ({ values }) => setConfig(normalizeConfig(specs, fromFormConfig(values))),
      }),
    [form, setConfig, specs],
  );

  const checks = check?.checks;
  const checked = check === undefined || checkPending ? undefined : settled;
  const configHash = runConfig?.status === "valid" ? runConfig.config_hash : undefined;

  const duplicate = useMemo(
    () =>
      configHash === undefined
        ? undefined
        : runList?.runs.find((run) => run.status === "ok" && run.command === command && run.config_hash === configHash),
    [command, configHash, runList],
  );

  const submit = useMemo(
    () => form.handleSubmit((values) => startRun(normalizeConfig(specs, fromFormConfig(values)))),
    [form, specs, startRun],
  );

  const onSubmit = useCallback((event: React.SubmitEvent) => void submit(event), [submit]);
  const run = useCallback(() => void submit(), [submit]);

  useHotkey(RUN_HOTKEY, run);

  const workspace = useMemo(
    (): WorkspaceView => ({ command, shownRunId, readOnly, settled, facts, check, run }),
    [check, command, facts, readOnly, run, settled, shownRunId],
  );

  return (
    <FormProvider {...form}>
      <WorkspaceContext value={workspace}>
        <WorkspaceLayout
          title={<WorkspaceTitle />}
          description={<WorkspaceDescription shownRunId={shownRunId} />}
          actions={<ResetButton />}
          panel={
            <form noValidate onSubmit={onSubmit} aria-label="Options">
              <OptionsPanel
                command={command}
                config={config}
                checks={checks}
                reported={check?.choices}
                checked={checked}
                facts={facts}
                effective={check?.effective}
              />
            </form>
          }
          footer={<RunFooter checks={checks} duplicate={duplicate} verb={COMMAND_VERBS[command]} onRun={run} />}
          main={<Outlet />}
        />
      </WorkspaceContext>
    </FormProvider>
  );
}

function WorkspaceTitle() {
  const label = useDraftStore((state) => datasetLabel(state.draft));

  return label ?? "New analysis";
}

function WorkspaceDescription({ shownRunId }: { shownRunId: JobId | undefined }) {
  const viewing = useDraftStore((state) => state.viewing?.runId);
  const fromRunId = useDraftStore((state) => state.draft.from_run_id);
  const { data: runList } = useApi((context) => runsList(context));
  const titleOf = (id: JobId) => runList?.runs.find((run) => run.id === id)?.title ?? id;

  const text =
    viewing === undefined
      ? fromRunId === undefined
        ? undefined
        : fromRunId === shownRunId
          ? "Variant of this run"
          : `Variant of "${titleOf(fromRunId)}"`
      : `Options of "${titleOf(viewing)}". An edit starts a variant of this run.`;

  return text === undefined ? null : <p className="text-muted-foreground text-xs">{text}</p>;
}

function ResetButton() {
  const command = useDraftStore((state) => state.draft.command);
  const discard = useDraftStore((state) => state.discard);
  const reset = useCallback(() => discard("reset", command), [command, discard]);

  return (
    <Button type="button" variant="ghost" size="sm" onClick={reset}>
      <RotateCcw aria-hidden />
      Reset form
    </Button>
  );
}
