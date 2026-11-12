import type { AppCommand, CheckConfigResponse, JobId } from "@neherlab/app-contracts";
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
import { COMMAND_LINE_HOTKEY, RUN_HOTKEY } from "../hotkeys";
import { commandSettings } from "../settings/catalog";
import { COMMAND_VERBS } from "../settings/commands";
import { normalizeConfig } from "../settings/config";
import { inputFactsRequest } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { RunFooter } from "./ChecksPanel";
import { type CommandLineCode, CommandLineDialog } from "./CommandLineDialog";
import { configResolver } from "./configResolver";
import { datasetLabel } from "./datasetLabel";
import { duplicateRun } from "./duplicateRun";
import { fromFormConfig, toFormConfig, type FormConfig } from "./formValues";
import { OptionsMenu } from "./OptionsMenu";
import { OptionsPanel } from "./OptionsPanel";
import { StaleBar } from "./StaleBar";
import { useStartRun } from "./useStartRun";
import { variantSummary } from "./variantSummary";
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
  const fromRunId = useDraftStore((state) => state.draft.from_run_id);
  const comparedRunId = shownRunId ?? fromRunId;

  const { data: check, isPlaceholderData: checkPending } = useApi(
    (context) => {
      const text = JSON.stringify(settled);
      const picked = choices ?? [];

      const body = { command, text, choices: picked };
      const withRun = comparedRunId === undefined ? body : { ...body, run_id: comparedRunId };

      return configCheck({
        ...context,
        body: facts === undefined ? withRun : { ...withRun, input_facts: facts },
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

  const duplicate = useMemo(
    () => duplicateRun(runList?.runs ?? [], command, runConfig?.status === "valid" ? runConfig : undefined),
    [command, runConfig, runList],
  );

  const submit = useMemo(
    () => form.handleSubmit((values) => startRun(normalizeConfig(specs, fromFormConfig(values)))),
    [form, specs, startRun],
  );

  const onSubmit = useCallback((event: React.SubmitEvent) => void submit(event), [submit]);
  const run = useCallback(() => void submit(), [submit]);

  useHotkey(RUN_HOTKEY, run);

  const [commandLineOpen, setCommandLineOpen] = useState(false);
  const showCommandLine = useCallback(() => setCommandLineOpen(true), []);

  useHotkey(COMMAND_LINE_HOTKEY, showCommandLine);

  const code = useMemo((): CommandLineCode => {
    if (runConfig === undefined) {
      return { state: "loading" };
    }

    return runConfig.status === "valid"
      ? { state: "ready", code: runConfig.code }
      : { state: "invalid", message: [runConfig.message, ...runConfig.causes].join(": ") };
  }, [runConfig]);

  const workspace = useMemo(
    (): WorkspaceView => ({ command, shownRunId, readOnly, settled, facts, check, run }),
    [check, command, facts, readOnly, run, settled, shownRunId],
  );

  return (
    <FormProvider {...form}>
      <WorkspaceContext value={workspace}>
        <WorkspaceLayout
          title={<WorkspaceTitle />}
          description={<WorkspaceDescription shownRunId={shownRunId} check={check} />}
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
                menu={<OptionsMenu command={command} onCommandLine={showCommandLine} />}
              />
            </form>
          }
          footer={<RunFooter checks={checks} duplicate={duplicate} verb={COMMAND_VERBS[command]} onRun={run} />}
          aboveMain={shownRunId === undefined ? undefined : <StaleBar shownRunId={shownRunId} />}
          main={<Outlet />}
        />
        <CommandLineDialog
          open={commandLineOpen}
          onOpenChange={setCommandLineOpen}
          description="The options as a command, with the paths that a run records"
          code={code}
        />
      </WorkspaceContext>
    </FormProvider>
  );
}

function WorkspaceTitle() {
  const label = useDraftStore((state) => datasetLabel(state.draft));

  return label ?? "New analysis";
}

function WorkspaceDescription({
  shownRunId,
  check,
}: {
  shownRunId: JobId | undefined;
  check: CheckConfigResponse | undefined;
}) {
  const viewing = useDraftStore((state) => state.viewing?.runId);
  const draft = useDraftStore((state) => state.draft);
  const { data: runList } = useApi((context) => runsList(context));
  const findRun = (id: JobId | undefined) => runList?.runs.find((run) => run.id === id);

  const text =
    viewing === undefined
      ? variantSummary({
          fromRunId: draft.from_run_id,
          shownRunId,
          viewing: false,
          source: findRun(draft.from_run_id),
          check,
          specs: commandSettings(draft.command).settings,
        })
      : `Options of "${findRun(viewing)?.title ?? viewing}". An edit starts a variant of this run.`;

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
