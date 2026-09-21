import * as z from "zod";

import type {
  AncestralArgs,
  AncestralResult,
  ClockArgs,
  ClockResult,
  DatasetInfo,
  LogEvent,
  MugrationArgs,
  MugrationResult,
  OptimizeArgs,
  OptimizeResult,
  ProgressEvent,
  PruneArgs,
  PruneResult,
  TimetreeArgs,
  TimetreeResult,
  VersionInfo,
} from "./generated/types.gen";
import {
  zAncestralResult,
  zClockResult,
  zDatasetInfo,
  zLogEvent,
  zMugrationResult,
  zOptimizeResult,
  zProgressEvent,
  zPruneResult,
  zTimetreeResult,
  zVersionInfo,
} from "./generated/zod.gen";

const zBridgeEvent = z.discriminatedUnion("type", [
  z.object({ type: z.literal("progress"), data: zProgressEvent }),
  z.object({ type: z.literal("log"), data: zLogEvent }),
]);

export type BridgeEvent = z.infer<typeof zBridgeEvent>;

export interface CommandOptions {
  onProgress?: (event: ProgressEvent) => void;
  signal?: AbortSignal;
}

export class CancelledError extends Error {
  constructor() {
    super("Operation cancelled");
    this.name = "CancelledError";
  }
}

export interface BridgeTransport {
  query(endpoint: string): Promise<unknown>;
  command(endpoint: string, args: unknown, options?: CommandOptions): Promise<unknown>;
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
    async version() {
      return zVersionInfo.parse(await transport.query("version"));
    },
    async datasets() {
      return z.array(zDatasetInfo).parse(await transport.query("datasets"));
    },
    async ancestral(args, options) {
      return zAncestralResult.parse(await transport.command("ancestral", args, options));
    },
    async clock(args, options) {
      return zClockResult.parse(await transport.command("clock", args, options));
    },
    async timetree(args, options) {
      return zTimetreeResult.parse(await transport.command("timetree", args, options));
    },
    async mugration(args, options) {
      return zMugrationResult.parse(await transport.command("mugration", args, options));
    },
    async optimize(args, options) {
      return zOptimizeResult.parse(await transport.command("optimize", args, options));
    },
    async prune(args, options) {
      return zPruneResult.parse(await transport.command("prune", args, options));
    },
  };
}

export function parseProgressEvent(data: unknown): ProgressEvent {
  return zProgressEvent.parse(data);
}

export function parseLogEvent(data: unknown): LogEvent {
  return zLogEvent.parse(data);
}

export function parseBridgeEvent(data: unknown): BridgeEvent {
  return zBridgeEvent.parse(data);
}
