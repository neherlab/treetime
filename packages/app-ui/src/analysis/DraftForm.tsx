import type { AppCommand } from "@neherlab/app-contracts";
import { configCheck, inputsCheck, runConfig as resolveRunConfig, runsList } from "@neherlab/app-contracts/client";
import { useHotkey } from "@tanstack/react-hotkeys";
import { useDebouncedValue } from "@tanstack/react-pacer";
import { keepPreviousData } from "@tanstack/react-query";
import { RotateCcw } from "lucide-react";
import { useCallback, useEffect, useMemo, useState } from "react";
import { FormProvider, useForm, useWatch } from "react-hook-form";

import { useApi } from "../api/hooks";
import { PageHeading, PageShell } from "../components/PageShell";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { changedSpecs, normalizeConfig } from "../settings/config";
import { inputFactsRequest } from "../settings/inputs";
import { zJsonObject } from "../settings/json";
import { useDraftStore } from "../store/draft";
import { Button } from "../ui/button";
import { ChecksPanel } from "./ChecksPanel";
import { CodePanel } from "./CodePanel";
import { CommandCards } from "./CommandCards";
import { configResolver } from "./configResolver";
import { toFormConfig, type FormConfig } from "./formValues";
import { InputsSection } from "./InputsSection";
import { SettingsPanel } from "./SettingsPanel";
import { useStartRun } from "./useStartRun";

const CHECK_DEBOUNCE = { wait: 400 };

export function DraftForm({ command }: { command: AppCommand }) {
  const specs = COMMAND_SETTINGS[command].specs;
  const [initial] = useState(() => normalizeConfig(specs, useDraftStore.getState().config));

  const form = useForm<FormConfig>({
    defaultValues: toFormConfig(initial),
    resolver: configResolver(command),
    mode: "onChange",
  });

  const setConfig = useDraftStore((state) => state.setConfig);
  const watched = useWatch({ control: form.control });

  const config = useMemo(() => {
    const parsed = zJsonObject.safeParse(watched);

    return parsed.success ? normalizeConfig(specs, parsed.data) : initial;
  }, [initial, specs, watched]);

  const [settled] = useDebouncedValue(config, CHECK_DEBOUNCE);
  const factsRequest = useMemo(() => inputFactsRequest(command, settled), [command, settled]);

  const { data: facts } = useApi(
    (context) => inputsCheck({ ...context, body: factsRequest ?? { command, config: {} } }),
    { enabled: factsRequest !== null, placeholderData: keepPreviousData, staleTime: 30_000 },
  );

  const { data: check } = useApi(
    (context) =>
      configCheck({ ...context, body: { command, text: JSON.stringify(settled), input_facts: facts ?? null } }),
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
        callback: ({ values }) => {
          const parsed = zJsonObject.safeParse(values);

          if (parsed.success) {
            setConfig(parsed.data);
          }
        },
      }),
    [form, setConfig],
  );

  const checks = check?.checks;
  const code = runConfig?.status === "valid" ? runConfig.code : null;

  const configHash = runConfig?.status === "valid" ? (runConfig.config_hash ?? null) : null;

  const duplicate = useMemo(
    () =>
      configHash === null
        ? undefined
        : runList?.runs.find((run) => run.status === "ok" && run.command === command && run.config_hash === configHash),
    [command, configHash, runList],
  );

  const changed = changedSpecs(specs, config);

  const submit = useMemo(
    () => form.handleSubmit((values) => startRun(normalizeConfig(specs, values))),
    [form, specs, startRun],
  );

  const onSubmit = useCallback((event: React.SubmitEvent) => void submit(event), [submit]);

  useHotkey("Mod+Enter", () => void submit());

  return (
    <PageShell className="max-w-[92.5rem]">
      <PageHeader />
      <FormProvider {...form}>
        <form
          noValidate
          onSubmit={onSubmit}
          className="grid grid-cols-[minmax(0,1fr)] items-start gap-5 @4xl:grid-cols-[minmax(0,1fr)_25rem]"
        >
          <div className="grid grid-cols-[minmax(0,1fr)] gap-6">
            <Step number={1} title="Analysis">
              <CommandCards command={command} config={config} />
            </Step>
            <Step number={2} title="Data">
              <InputsSection command={command} facts={facts ?? undefined} />
            </Step>
            <Step
              number={3}
              title="Settings"
              hint={changed.length === 0 ? "All defaults" : `${changed.length} changed`}
            >
              <SettingsPanel command={command} config={config} facts={facts ?? undefined} checks={checks} />
            </Step>
          </div>
          <aside className="grid gap-4 @4xl:sticky @4xl:top-5">
            <CodePanel command={command} code={code} />
            <ChecksPanel checks={checks} duplicate={duplicate} verb={COMMAND_INFO[command].verb} />
          </aside>
        </form>
      </FormProvider>
    </PageShell>
  );
}

function PageHeader() {
  const fromRunId = useDraftStore((state) => state.fromRunId);
  const command = useDraftStore((state) => state.command);
  const reset = useDraftStore((state) => state.reset);
  const { data: runList } = useApi((context) => runsList(context));
  const fromTitle = runList?.runs.find((run) => run.id === fromRunId)?.title ?? fromRunId;
  const resetForm = useCallback(() => reset(command), [command, reset]);

  return (
    <PageHeading
      title={fromRunId === null ? "New analysis" : "Edit and run again"}
      description={
        fromRunId === null
          ? "Choose an analysis, add the data, check the settings, run."
          : `Settings copied from "${fromTitle ?? ""}". The original run stays unchanged.`
      }
      actions={
        <Button type="button" variant="ghost" onClick={resetForm}>
          <RotateCcw aria-hidden />
          Reset form
        </Button>
      }
    />
  );
}

function Step({
  number,
  title,
  hint,
  children,
}: {
  number: number;
  title: string;
  hint?: string;
  children: React.ReactNode;
}) {
  return (
    <section aria-labelledby={`step-${number}`} className="@container grid gap-2.5">
      <div className="flex items-baseline gap-2.5 border-b pb-1.5">
        <span className="text-primary font-mono text-sm font-medium" aria-hidden>
          {String(number).padStart(2, "0")}
        </span>
        <h2 id={`step-${number}`} className="font-heading text-base font-semibold">
          {title}
        </h2>
        {hint !== undefined && <span className="text-muted-foreground text-sm">{hint}</span>}
      </div>
      {children}
    </section>
  );
}
