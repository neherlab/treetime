import type {
  AncestralArgs,
  AncestralResult,
  ClockArgs,
  ClockResult,
  MugrationArgs,
  MugrationResult,
  OptimizeArgs,
  OptimizeResult,
  PruneArgs,
  PruneResult,
  TimetreeArgs,
  TimetreeResult,
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

export function useAncestral() {
  const bridge = useBridge();
  return useMutation<AncestralResult, Error, AncestralArgs>({
    mutationFn: (args) => bridge.ancestral(args),
  });
}

export function useClock() {
  const bridge = useBridge();
  return useMutation<ClockResult, Error, ClockArgs>({
    mutationFn: (args) => bridge.clock(args),
  });
}

export function useTimetree() {
  const bridge = useBridge();
  return useMutation<TimetreeResult, Error, TimetreeArgs>({
    mutationFn: (args) => bridge.timetree(args),
  });
}

export function useMugration() {
  const bridge = useBridge();
  return useMutation<MugrationResult, Error, MugrationArgs>({
    mutationFn: (args) => bridge.mugration(args),
  });
}

export function useOptimize() {
  const bridge = useBridge();
  return useMutation<OptimizeResult, Error, OptimizeArgs>({
    mutationFn: (args) => bridge.optimize(args),
  });
}

export function usePrune() {
  const bridge = useBridge();
  return useMutation<PruneResult, Error, PruneArgs>({
    mutationFn: (args) => bridge.prune(args),
  });
}
