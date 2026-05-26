import type { AncestralArgs, ClockArgs, MugrationArgs, OptimizeArgs, PruneArgs, TimetreeArgs } from "./args";
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
import type { CommandResult, ProgressEvent } from "./results";
=======
import type { DatasetInfo } from "./datasets";
import type {
  AncestralResult,
  ClockResult,
=======
=======
import type { DatasetInfo } from "./datasets";
>>>>>>> b8625b9a (feat(app): wire datasets through bridge contract and implementations)
import type {
  AncestralResult,
  ClockResult,
  HomoplasyResult,
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
  MugrationResult,
  OptimizeResult,
  ProgressEvent,
  PruneResult,
  TimetreeResult,
} from "./results";
<<<<<<< HEAD
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)
=======
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
import type { VersionInfo } from "./version";

export interface CommandOptions {
  onProgress?: (event: ProgressEvent) => void;
  signal?: AbortSignal;
}

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
  ancestral(args: AncestralArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  clock(args: ClockArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  timetree(args: TimetreeArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  mugration(args: MugrationArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  optimize(args: OptimizeArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  prune(args: PruneArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
=======
  datasets(): Promise<DatasetInfo[]>;
  ancestral(args: AncestralArgs, options?: CommandOptions): Promise<AncestralResult>;
  clock(args: ClockArgs, options?: CommandOptions): Promise<ClockResult>;
  timetree(args: TimetreeArgs, options?: CommandOptions): Promise<TimetreeResult>;
  mugration(args: MugrationArgs, options?: CommandOptions): Promise<MugrationResult>;
  optimize(args: OptimizeArgs, options?: CommandOptions): Promise<OptimizeResult>;
  prune(args: PruneArgs, options?: CommandOptions): Promise<PruneResult>;
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)
=======
=======
  datasets(): Promise<DatasetInfo[]>;
>>>>>>> b8625b9a (feat(app): wire datasets through bridge contract and implementations)
  ancestral(args: AncestralArgs, onProgress?: (event: ProgressEvent) => void): Promise<AncestralResult>;
  clock(args: ClockArgs, onProgress?: (event: ProgressEvent) => void): Promise<ClockResult>;
  timetree(args: TimetreeArgs, onProgress?: (event: ProgressEvent) => void): Promise<TimetreeResult>;
  mugration(args: MugrationArgs, onProgress?: (event: ProgressEvent) => void): Promise<MugrationResult>;
  optimize(args: OptimizeArgs, onProgress?: (event: ProgressEvent) => void): Promise<OptimizeResult>;
  prune(args: PruneArgs, onProgress?: (event: ProgressEvent) => void): Promise<PruneResult>;
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
}
