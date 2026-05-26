// Generated from openapi.yaml - do not edit

import type { AncestralArgs, AncestralResult, ClockArgs, ClockResult, DatasetInfo, MugrationArgs, MugrationResult, OptimizeArgs, OptimizeResult, PruneArgs, PruneResult, TimetreeArgs, TimetreeResult, VersionInfo } from "./types";

import type { ProgressEvent } from "./types";

export interface CommandOptions {
  onProgress?: (event: ProgressEvent) => void;
  signal?: AbortSignal;
}

export class CancelledError extends Error {
  constructor() { super("Operation cancelled"); this.name = "CancelledError"; }
}

export interface BridgeTransport {
  query<T>(endpoint: string): Promise<T>;
  command<T>(endpoint: string, args: unknown, options?: CommandOptions): Promise<T>;
}

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
  datasets(): Promise<DatasetInfo[]>;
  ancestral(args: AncestralArgs, options?: CommandOptions): Promise<AncestralResult>;
  clock(args: ClockArgs, options?: CommandOptions): Promise<ClockResult>;
  timetree(args: TimetreeArgs, options?: CommandOptions): Promise<TimetreeResult>;
  mugration(args: MugrationArgs, options?: CommandOptions): Promise<MugrationResult>;
  optimize(args: OptimizeArgs, options?: CommandOptions): Promise<OptimizeResult>;
  prune(args: PruneArgs, options?: CommandOptions): Promise<PruneResult>;
}

export function createBridge(transport: BridgeTransport): TreeTimeBridge {
  return {
    version: () => transport.query<VersionInfo>("version"),
    datasets: () => transport.query<DatasetInfo[]>("datasets"),
    ancestral: (args, options) => transport.command<AncestralResult>("ancestral", args, options),
    clock: (args, options) => transport.command<ClockResult>("clock", args, options),
    timetree: (args, options) => transport.command<TimetreeResult>("timetree", args, options),
    mugration: (args, options) => transport.command<MugrationResult>("mugration", args, options),
    optimize: (args, options) => transport.command<OptimizeResult>("optimize", args, options),
    prune: (args, options) => transport.command<PruneResult>("prune", args, options),
  };
}
