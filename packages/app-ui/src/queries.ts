import type { AppCommand, CheckInputsRequest } from "@neherlab/app-contracts";
import { keepPreviousData, useQueries, useQuery } from "@tanstack/react-query";

import { useBridge } from "./BridgeContext";
import { canonicalJson, type JsonObject } from "./settings/json";

const ACTIVE_RUNS_POLL_MS = 2000;

const IDLE_RUNS_POLL_MS = 15_000;

export const RUNS_KEY = ["runs"] as const;

export function useVersion() {
  const bridge = useBridge();

  return useQuery({ queryKey: ["version"], queryFn: () => bridge.version(), staleTime: Infinity });
}

export function useDatasetCatalog() {
  const bridge = useBridge();

  return useQuery({ queryKey: ["datasets"], queryFn: () => bridge.datasets(), staleTime: Infinity });
}

export function useRunList() {
  const bridge = useBridge();

  return useQuery({
    queryKey: RUNS_KEY,
    queryFn: () => bridge.listRuns(),
    staleTime: 0,
    refetchInterval: (query) => ((query.state.data?.active_runs ?? 0) > 0 ? ACTIVE_RUNS_POLL_MS : IDLE_RUNS_POLL_MS),
  });
}

export function useRunRecord(id: string) {
  const bridge = useBridge();

  return useQuery({
    queryKey: [...RUNS_KEY, id],
    queryFn: () => bridge.getRun(id),
    staleTime: 0,
    refetchInterval: (query) => (query.state.data?.status === "running" ? ACTIVE_RUNS_POLL_MS : false),
  });
}

export function useRunRecords(ids: readonly string[]) {
  const bridge = useBridge();

  return useQueries({
    queries: ids.map((id) => ({ queryKey: [...RUNS_KEY, id], queryFn: () => bridge.getRun(id), staleTime: 60_000 })),
  });
}

export function useConfigCheck(command: AppCommand, config: JsonObject) {
  const bridge = useBridge();
  const text = JSON.stringify(config);

  return useQuery({
    queryKey: ["check-config", command, text],
    queryFn: () => bridge.checkConfig({ command, text }),
    placeholderData: keepPreviousData,
    staleTime: Infinity,
  });
}

export function useRunConfig(command: AppCommand, config: JsonObject) {
  const bridge = useBridge();

  return useQuery({
    queryKey: ["run-config", command, canonicalJson(config)],
    queryFn: () => bridge.runConfig({ command, config }),
    placeholderData: keepPreviousData,
    staleTime: 10_000,
  });
}

export function useInputFacts(request: CheckInputsRequest | null) {
  const bridge = useBridge();

  return useQuery({
    queryKey: ["check-inputs", request],
    queryFn: () => (request === null ? Promise.resolve(null) : bridge.checkInputs(request)),
    enabled: request !== null,
    placeholderData: keepPreviousData,
    staleTime: 30_000,
  });
}
