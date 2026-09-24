import * as z from "zod";

import type {
  AncestralConfig,
  AppCommand,
  CheckConfigRequest,
  ClockConfig,
  CommandOutcome,
  DatasetInfo,
  JobEvent,
  LogEvent,
  MugrationConfig,
  OptimizeConfig,
  ProgressEvent,
  PruneConfig,
  TerminalEvent,
  TimetreeConfig,
  VersionInfo,
} from "./generated/types.gen";
import { zCheckConfigResponse, zDatasetInfo, zJobEvent, zTerminalEvent, zVersionInfo } from "./generated/zod.gen";

export type CheckConfigResult = z.infer<typeof zCheckConfigResponse>;

export interface CommandOptions {
  onStarted?: (jobId: string) => void;
  onProgress?: (event: ProgressEvent) => void;
  onLog?: (event: LogEvent) => void;
  signal?: AbortSignal;
}

export interface TransportCommandOptions {
  onEvent: (event: JobEvent) => void;
  signal?: AbortSignal;
}

export class CancelledError extends Error {
  constructor() {
    super("Operation cancelled");
    this.name = "CancelledError";
  }
}

export class CommandError extends Error {
  readonly jobId: string;
  readonly causes: string[];

  constructor(jobId: string, message: string, causes: string[]) {
    super(message);
    this.name = "CommandError";
    this.jobId = jobId;
    this.causes = causes;
  }
}

export interface BridgeTransport {
  query(endpoint: string): Promise<unknown>;
  request(endpoint: string, body: unknown): Promise<unknown>;
  command(command: AppCommand, config: unknown, options: TransportCommandOptions): Promise<unknown>;
}

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
  datasets(): Promise<DatasetInfo[]>;
  checkConfig(request: CheckConfigRequest): Promise<CheckConfigResult>;
  timetree(config: TimetreeConfig, options?: CommandOptions): Promise<CommandOutcome>;
  optimize(config: OptimizeConfig, options?: CommandOptions): Promise<CommandOutcome>;
  prune(config: PruneConfig, options?: CommandOptions): Promise<CommandOutcome>;
  ancestral(config: AncestralConfig, options?: CommandOptions): Promise<CommandOutcome>;
  clock(config: ClockConfig, options?: CommandOptions): Promise<CommandOutcome>;
  mugration(config: MugrationConfig, options?: CommandOptions): Promise<CommandOutcome>;
}

export function createBridge(transport: BridgeTransport): TreeTimeBridge {
  async function run(command: AppCommand, config: unknown, options: CommandOptions = {}): Promise<CommandOutcome> {
    const onEvent = (event: JobEvent) => {
      dispatchEvent(event, options);
    };

    const transportOptions: TransportCommandOptions = { onEvent };

    if (options.signal !== undefined) {
      transportOptions.signal = options.signal;
    }

    const terminal = parseTerminalEvent(await transport.command(command, config, transportOptions));

    return commandOutcome(terminal);
  }

  return {
    async version() {
      return zVersionInfo.parse(await transport.query("version"));
    },
    async datasets() {
      return z.array(zDatasetInfo).parse(await transport.query("datasets"));
    },
    async checkConfig(request) {
      return zCheckConfigResponse.parse(await transport.request("check-config", request));
    },
    timetree: (config, options) => run("timetree", config, options),
    optimize: (config, options) => run("optimize", config, options),
    prune: (config, options) => run("prune", config, options),
    ancestral: (config, options) => run("ancestral", config, options),
    clock: (config, options) => run("clock", config, options),
    mugration: (config, options) => run("mugration", config, options),
  };
}

function commandOutcome(terminal: TerminalEvent): CommandOutcome {
  if (terminal.status === "error") {
    throw new CommandError(terminal.job_id, terminal.message, terminal.causes);
  }

  if (terminal.status === "cancelled") {
    throw new CancelledError();
  }

  return terminal.result;
}

function dispatchEvent(event: JobEvent, options: CommandOptions): void {
  switch (event.type) {
    case "started":
      options.onStarted?.(event.data.job_id);
      break;
    case "progress":
      options.onProgress?.(event.data);
      break;
    case "log":
      options.onLog?.(event.data);
      break;
    case "terminal":
      break;
  }
}

export function parseJobEvent(data: unknown): JobEvent {
  return zJobEvent.parse(data);
}

function parseTerminalEvent(data: unknown): TerminalEvent {
  return zTerminalEvent.parse(data);
}
