import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback, useEffect, useState } from "react";
import { type Resolver, useForm } from "react-hook-form";

import { commandSettings } from "../settings/catalog";
import { normalizeConfig } from "../settings/config";
import { useDraftStore } from "../store/draft";
import { configResolver } from "./configResolver";
import { fromFormConfig, toFormConfig, type FormConfig } from "./formValues";

export function useDraftForm() {
  const [defaultValues] = useState(draftFormValues);

  const resolver = useCallback<Resolver<FormConfig>>(
    (values, context, options) => configResolver(currentCommand())(values, context, options),
    [],
  );

  const form = useForm<FormConfig>({ defaultValues, resolver, mode: "onChange" });

  useEffect(
    () =>
      useDraftStore.subscribe((state, previous) => {
        if (state.epoch !== previous.epoch) {
          form.reset(draftFormValues());
        }
      }),
    [form],
  );

  useEffect(
    () =>
      form.subscribe({
        formState: { values: true },
        callback: ({ values }) => {
          const specs = commandSettings(currentCommand()).settings;

          useDraftStore.getState().setConfig(normalizeConfig(specs, fromFormConfig(values)));
        },
      }),
    [form],
  );

  return form;
}

function draftFormValues(): FormConfig {
  const { draft } = useDraftStore.getState();

  return toFormConfig(normalizeConfig(commandSettings(draft.command).settings, draft.config));
}

function currentCommand(): AppCommand {
  return useDraftStore.getState().draft.command;
}
