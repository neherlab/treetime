import type {
  AncestralConfig,
  ClockConfig,
  CommandOutcome,
  DatasetInfo,
  MugrationConfig,
  OptimizeConfig,
  PruneConfig,
  TimetreeConfig,
  VersionInfo,
} from "@neherlab/app-contracts";
import { useMutation, useQuery } from "@tanstack/react-query";

import { useBridge } from "./BridgeContext";

export function useVersion() {
  const bridge = useBridge();

  return useQuery<VersionInfo>({
    queryKey: ["version"],
    queryFn: () => bridge.version(),
    staleTime: Infinity,
  });
}

export function useDatasets() {
  const bridge = useBridge();

  return useQuery<DatasetInfo[]>({
    queryKey: ["datasets"],
    queryFn: () => bridge.datasets(),
    staleTime: Infinity,
  });
}

export function useAncestral() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, AncestralConfig>({
    mutationFn: (args) => bridge.ancestral(args),
  });
}

export function useClock() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, ClockConfig>({
    mutationFn: (args) => bridge.clock(args),
  });
}

export function useTimetree() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, TimetreeConfig>({
    mutationFn: (args) => bridge.timetree(args),
  });
}

export function useMugration() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, MugrationConfig>({
    mutationFn: (args) => bridge.mugration(args),
  });
}

export function useOptimize() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, OptimizeConfig>({
    mutationFn: (args) => bridge.optimize(args),
  });
}

export function usePrune() {
  const bridge = useBridge();

  return useMutation<CommandOutcome, Error, PruneConfig>({
    mutationFn: (args) => bridge.prune(args),
  });
}
