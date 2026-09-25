import type { AppCommand, CheckInputsRequest, InputFactsResult, TreeTimeBridge } from "@neherlab/app-contracts";
import { keepPreviousData, useQueries, useQuery } from "@tanstack/react-query";

import { useBridge } from "./BridgeContext";
import { indexClades } from "./results/clades";
import { outputPath } from "./results/files";
import { loadRunResults } from "./results/load";
import { parseAuspiceJson, readAuspiceTree } from "./results/tree";
import { canonicalJson, type JsonObject } from "./settings/json";

const ACTIVE_RUNS_POLL_MS = 2000;

const IDLE_RUNS_POLL_MS = 15_000;

export const RUNS_KEY = ["runs"] as const;

const OUTPUTS_KEY = "run-outputs";

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

export function useConfigCheck(command: AppCommand, config: JsonObject, facts: InputFactsResult | null | undefined) {
  const bridge = useBridge();
  const text = JSON.stringify(config);
  const inputFacts = facts ?? null;

  return useQuery({
    queryKey: ["check-config", command, text, inputFacts],
    queryFn: () => bridge.checkConfig({ command, text, input_facts: inputFacts }),
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

export function useRunFiles(id: string, enabled: boolean) {
  const bridge = useBridge();

  return useQuery({
    queryKey: [OUTPUTS_KEY, id, "files"],
    queryFn: () => bridge.runFiles(id),
    enabled,
    staleTime: Infinity,
  });
}

export function useRunResults(id: string, command: AppCommand, enabled: boolean) {
  const bridge = useBridge();

  return useQuery({
    queryKey: [OUTPUTS_KEY, id, "results"],
    queryFn: async () => loadRunResults(command, await bridge.runFiles(id), (path) => readText(bridge, id, path)),
    enabled,
    staleTime: Infinity,
  });
}

export function useRunTrees(ids: readonly string[]) {
  const bridge = useBridge();

  return useQueries({
    queries: ids.map((id) => ({
      queryKey: [OUTPUTS_KEY, id, "tree"],
      queryFn: async () => {
        const path = outputPath(await bridge.runFiles(id), "auspice");

        if (path === undefined) {
          return null;
        }

        const tree = readAuspiceTree(parseAuspiceJson(await readText(bridge, id, path)));

        return { tree, clades: indexClades(tree) };
      },
      staleTime: Infinity,
    })),
  });
}

async function readText(bridge: TreeTimeBridge, id: string, path: string): Promise<string> {
  return new TextDecoder().decode(await bridge.readRunFile(id, path));
}
