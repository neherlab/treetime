import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback, useEffect, useMemo, useState } from "react";
import { FormProvider, useForm, useWatch } from "react-hook-form";
import { useDebounce } from "use-debounce";

import { useConfigCheck, useInputFacts, useRunConfig, useRunList } from "../queries";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { changedSpecs, normalizeConfig, outputFreeConfig } from "../settings/config";
import { inputFactsRequest } from "../settings/inputs";
import { zJsonObject } from "../settings/json";
import { autoTitle, changedChips } from "../settings/titles";
import { useDraftStore } from "../store/draft";
import { ChecksPanel } from "./ChecksPanel";
import { CodePanel } from "./CodePanel";
import { CommandCards } from "./CommandCards";
import { configResolver } from "./configResolver";
import { toFormConfig, type FormConfig } from "./formValues";
import { InputsSection } from "./InputsSection";
import { SettingsPanel } from "./SettingsPanel";
import { useStartRun } from "./useStartRun";

const CHECK_DEBOUNCE_MS = 400;

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

  const [settled] = useDebounce(config, CHECK_DEBOUNCE_MS);
  const request = useMemo(() => outputFreeConfig(specs, settled), [settled, specs]);
  const { data: facts } = useInputFacts(useMemo(() => inputFactsRequest(command, settled), [command, settled]));
  const { data: configCheck } = useConfigCheck(command, request, facts);
  const { data: runConfig } = useRunConfig(command, request);
  const { data: runList } = useRunList();
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

  const checks = configCheck?.checks;
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

  const title = autoTitle(
    command,
    changedChips(
      code,
      specs,
      changed.map((spec) => spec.key),
    ),
  );

  const submit = form.handleSubmit((values) =>
    startRun(outputFreeConfig(specs, normalizeConfig(specs, values)), title),
  );

  const onSubmit = useCallback((event: React.FormEvent) => void submit(event), [submit]);

  useEffect(() => {
    function onKeyDown(event: KeyboardEvent) {
      if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
        event.preventDefault();
        void submit();
      }
    }

    window.addEventListener("keydown", onKeyDown);

    return () => window.removeEventListener("keydown", onKeyDown);
  }, [submit]);

  return (
    <div className="mx-auto max-w-[92.5rem] px-5 pt-4 pb-16">
      <PageHeader />
      <FormProvider {...form}>
        <form noValidate onSubmit={onSubmit} className="grid items-start gap-4 xl:grid-cols-[minmax(0,1fr)_25rem]">
          <div className="grid min-w-0 gap-5">
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
              <SettingsPanel command={command} config={config} facts={facts ?? undefined} />
            </Step>
          </div>
          <aside className="grid gap-3.5 xl:sticky xl:top-4">
            <CodePanel command={command} code={code} />
            <ChecksPanel
              suggestedTitle={title}
              checks={checks}
              duplicate={duplicate}
              verb={COMMAND_INFO[command].verb}
            />
          </aside>
        </form>
      </FormProvider>
    </div>
  );
}

function PageHeader() {
  const fromRunId = useDraftStore((state) => state.fromRunId);
  const command = useDraftStore((state) => state.command);
  const reset = useDraftStore((state) => state.reset);
  const { data: runList } = useRunList();
  const fromTitle = runList?.runs.find((run) => run.id === fromRunId)?.title ?? fromRunId;
  const resetForm = useCallback(() => reset(command), [command, reset]);

  return (
    <div className="mb-3.5 flex items-start gap-4">
      <div>
        <h1 className="text-2xl leading-tight font-bold">
          {fromRunId === null ? "New analysis" : "Edit and run again"}
        </h1>
        <p className="text-ink-muted mt-1">
          {fromRunId === null
            ? "Choose an analysis, add the data, check the settings, run."
            : `Settings copied from "${fromTitle ?? ""}". The original run stays unchanged.`}
        </p>
      </div>
      <button
        type="button"
        onClick={resetForm}
        className="text-ink-muted hover:bg-surface-2 hover:text-ink ml-auto rounded-md px-3 py-1.5 font-bold"
      >
        Reset form
      </button>
    </div>
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
    <section aria-labelledby={`step-${number}`}>
      <div className="mb-2 flex items-baseline gap-2.5">
        <span className="bg-ink text-surface-1 inline-grid size-5.5 place-items-center rounded-full text-xs font-bold">
          {number}
        </span>
        <h2 id={`step-${number}`} className="text-base font-bold">
          {title}
        </h2>
        {hint !== undefined && <span className="text-ink-faint">{hint}</span>}
      </div>
      {children}
    </section>
  );
}
