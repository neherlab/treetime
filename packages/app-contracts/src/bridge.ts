import type { AncestralArgs, ClockArgs, MugrationArgs, OptimizeArgs, PruneArgs, TimetreeArgs } from "./args";
import type { CommandResult, ProgressEvent } from "./results";
import type { VersionInfo } from "./version";

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
  ancestral(args: AncestralArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  clock(args: ClockArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  timetree(args: TimetreeArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  mugration(args: MugrationArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  optimize(args: OptimizeArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
  prune(args: PruneArgs, onProgress?: (event: ProgressEvent) => void): Promise<CommandResult>;
}
