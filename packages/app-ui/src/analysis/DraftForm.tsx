import type { AppCommand, JobId, RunCheck } from "@neherlab/app-contracts";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useHotkey } from "@tanstack/react-hotkeys";
import { Outlet } from "@tanstack/react-router";
import { useCallback, useMemo, useState } from "react";
import { FormProvider, useFormState, useWatch } from "react-hook-form";

import { LoadingState } from "../components/PageShell";
import { COMMAND_LINE_HOTKEY, RUN_HOTKEY } from "../hotkeys";
import { commandSettings } from "../settings/catalog";
import { COMMAND_VERBS } from "../settings/commands";
import { normalizeConfig } from "../settings/config";
import { useListedRuns } from "../shell/useListedRuns";
import { useDraftStore } from "../store/draft";
import { RunFooter } from "./ChecksPanel";
import { type CommandLineCode, CommandLineDialog } from "./CommandLineDialog";
import { duplicateRun } from "./duplicateRun";
import { fromFormConfig, type FormConfig } from "./formValues";
import { OptionsMenu } from "./OptionsMenu";
import { OptionsPanel } from "./OptionsPanel";
import { RUN_BLOCK_TEXT, runBlock, SingleRun } from "./runAction";
import { StaleBar } from "./StaleBar";
import { useDraftChecks } from "./useDraftChecks";
import { useDraftForm } from "./useDraftForm";
import { useStartRun } from "./useStartRun";
import { WorkspaceContext, type WorkspaceView } from "./workspaceContext";
import { ResetButton, WorkspaceDescription, WorkspaceTitle } from "./WorkspaceHeader";
import { WorkspaceLayout } from "./WorkspaceLayout";

export function DraftForm({
  command,
  shownRunId,
  pending,
}: {
  command: AppCommand;
  shownRunId: JobId | undefined;
  pending: boolean;
}) {
  const specs = commandSettings(command).settings;
  const form = useDraftForm();
  const watched = useWatch({ control: form.control });
  const config = useMemo(() => normalizeConfig(specs, fromFormConfig(watched)), [specs, watched]);
  const { settled, facts, check, checkPending, runConfig } = useDraftChecks(command, config, shownRunId);
  const { all: listedRuns } = useListedRuns();
  const run = useRunAction(form, command, check?.checks, checkPending);
  const epoch = useDraftStore((state) => state.epoch);

  const duplicate = useMemo(
    () => duplicateRun(listedRuns, command, runConfig?.status === "valid" ? runConfig : undefined),
    [command, listedRuns, runConfig],
  );

  const [commandLineOpen, setCommandLineOpen] = useState(false);
  const showCommandLine = useCallback(() => setCommandLineOpen(true), []);

  useHotkey(COMMAND_LINE_HOTKEY, showCommandLine);

  const workspace = useMemo(
    (): WorkspaceView => ({ command, settled, facts, check, checkPending, ...run }),
    [check, checkPending, command, facts, run, settled],
  );

  return (
    <FormProvider {...form}>
      <WorkspaceContext value={workspace}>
        <WorkspaceLayout
          title={<WorkspaceTitle check={check} />}
          description={<WorkspaceDescription shownRunId={shownRunId} check={check} />}
          actions={<ResetButton />}
          panel={
            pending ? (
              <LoadingState text="Loading the options of the run" />
            ) : (
              <OptionsPanel
                key={`${command}:${epoch}`}
                command={command}
                config={config}
                check={check}
                checked={check === undefined || checkPending ? undefined : settled}
                facts={facts}
                menu={<OptionsMenu command={command} onCommandLine={showCommandLine} />}
              />
            )
          }
          footer={<RunFooter checks={check?.checks} duplicate={duplicate} verb={COMMAND_VERBS[command]} />}
          aboveMain={shownRunId === undefined ? undefined : <StaleBar shownRunId={shownRunId} />}
          main={<Outlet />}
        />
        <CommandLineDialog
          open={commandLineOpen}
          onOpenChange={setCommandLineOpen}
          description="The options as a command, with the paths that a run records"
          code={commandLineCode(runConfig)}
        />
      </WorkspaceContext>
    </FormProvider>
  );
}

function useRunAction(
  form: ReturnType<typeof useDraftForm>,
  command: AppCommand,
  checks: readonly RunCheck[] | undefined,
  checkPending: boolean,
) {
  const startRun = useStartRun(command);
  const toasts = useToastManager();
  const { isValid } = useFormState<FormConfig>({ control: form.control });
  const [single] = useState(() => new SingleRun());
  const [starting, setStarting] = useState(false);
  const block = runBlock({ starting, checks, checkPending, valid: isValid });

  const start = useCallback(async () => {
    setStarting(true);

    const submit = form.handleSubmit(async (values) =>
      startRun(normalizeConfig(commandSettings(command).settings, fromFormConfig(values))),
    );

    await submit().finally(() => {
      setStarting(false);
    });
  }, [command, form, startRun]);

  const run = useCallback(() => {
    if (block !== undefined) {
      toasts.add({ title: "The run cannot start yet", description: RUN_BLOCK_TEXT[block] });

      return;
    }

    void single.run(start);
  }, [block, single, start, toasts]);

  useHotkey(RUN_HOTKEY, run);

  return useMemo(() => ({ run, runBlock: block }), [block, run]);
}

function commandLineCode(runConfig: ReturnType<typeof useDraftChecks>["runConfig"]): CommandLineCode {
  if (runConfig === undefined) {
    return { state: "loading" };
  }

  return runConfig.status === "valid"
    ? { state: "ready", code: runConfig.code }
    : { state: "invalid", message: [runConfig.message, ...runConfig.causes].join(": ") };
}
