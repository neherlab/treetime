import type { AppCommand, CheckConfigRequest, JobId, SparseConfig } from "@neherlab/app-contracts";
import { configCheck, inputsCheck, runConfig as resolveRunConfig } from "@neherlab/app-contracts/client";
import { useDebouncedValue } from "@tanstack/react-pacer";
import { keepPreviousData } from "@tanstack/react-query";
import { useMemo } from "react";

import { useApi } from "../api/hooks";
import { inputFactsRequest } from "../settings/inputs";
import { useDraftStore } from "../store/draft";

const CHECK_DEBOUNCE = { wait: 400 };

export function useDraftChecks(command: AppCommand, config: SparseConfig, shownRunId: JobId | undefined) {
  const [settled] = useDebouncedValue(config, CHECK_DEBOUNCE);
  const factsRequest = useMemo(() => inputFactsRequest(command, settled), [command, settled]);

  const { data: facts } = useApi(
    (context) => inputsCheck({ ...context, body: factsRequest ?? { command, config: {} } }),
    { enabled: factsRequest !== null, placeholderData: keepPreviousData, staleTime: 30_000 },
  );

  const choices = useDraftStore((state) => state.draft.choices);
  const fromRunId = useDraftStore((state) => state.draft.from_run_id);
  const example = useDraftStore((state) => state.draft.example);
  const comparedRunId = shownRunId ?? fromRunId;

  const { data: check, isPlaceholderData: checkPending } = useApi(
    (context) => {
      const body: CheckConfigRequest = { command, text: JSON.stringify(settled), choices: choices ?? [] };

      if (comparedRunId !== undefined) {
        body.run_id = comparedRunId;
      }

      if (example !== undefined) {
        body.example = example;
      }

      if (facts !== undefined) {
        body.input_facts = facts;
      }

      return configCheck({ ...context, body });
    },
    { placeholderData: keepPreviousData, staleTime: Infinity },
  );

  const { data: runConfig } = useApi(
    (context) => resolveRunConfig({ ...context, body: { command, config: settled } }),
    { placeholderData: keepPreviousData, staleTime: 10_000 },
  );

  return { settled, facts, check, checkPending, runConfig };
}
