import * as z from "zod";

import type {
  AncestralConfig,
  AppCommand,
  CheckConfigRequest,
  CheckInputsRequest,
  ClockConfig,
  CreateRunRequest,
  DatasetCatalog,
  LogEvent,
  MugrationConfig,
  OptimizeConfig,
  ProgressEvent,
  PruneConfig,
  StartRunRequest,
  TimetreeConfig,
  UpdateRunRequest,
  VersionInfo,
} from "./generated/types.gen";
import {
  zCheckConfigResponse,
  zCommandOutcome,
  zDatasetCatalog,
  zInputFacts,
  zIterationEvent,
  zRunEvent,
  zRunFile,
  zRunList,
  zRunRecord,
  zRunSummary,
  zTerminalEvent,
  zUploadedInput,
  zVersionInfo,
} from "./generated/zod.gen";

export type Parsed<S extends z.ZodType> = z.infer<S>;

export type CheckConfigResult = Parsed<typeof zCheckConfigResponse>;

type RunEventResult = Parsed<typeof zRunEvent>;

type TerminalEventResult = Parsed<typeof zTerminalEvent>;

type CommandOutcomeResult = Parsed<typeof zCommandOutcome>;

export interface CommandOptions {
  title?: string;
  onStarted?: (runId: string) => void;
  onProgress?: (event: ProgressEvent) => void;
  onLog?: (event: LogEvent) => void;
  onIteration?: (event: Parsed<typeof zIterationEvent>) => void;
  signal?: AbortSignal;
}

export interface FollowRunOptions {
  from?: number;
  onEvent?: (event: RunEventResult) => void;
  signal?: AbortSignal;
}

export interface TransportEventOptions {
  from: number;
  onEvent: (event: unknown) => void;
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

export class RunEndedError extends Error {
  readonly runId: string;

  constructor(runId: string) {
    super(`the event stream of run ${runId} ended without a terminal event`);
    this.name = "RunEndedError";
    this.runId = runId;
  }
}

export interface BridgeTransport {
  version(): Promise<unknown>;
  datasets(): Promise<unknown>;
  checkConfig(request: CheckConfigRequest): Promise<unknown>;
  checkInputs(request: CheckInputsRequest): Promise<unknown>;
  listRuns(): Promise<unknown>;
  createRun(request: CreateRunRequest): Promise<unknown>;
  getRun(id: string): Promise<unknown>;
  startRun(id: string, request: StartRunRequest): Promise<unknown>;
  updateRun(id: string, request: UpdateRunRequest): Promise<unknown>;
  cancelRun(id: string): Promise<unknown>;
  deleteRun(id: string): Promise<void>;
  restoreRun(id: string): Promise<unknown>;
  purgeRun(id: string): Promise<void>;
  runEvents(id: string, options: TransportEventOptions): Promise<void>;
  runFiles(id: string): Promise<unknown>;
  readRunFile(id: string, path: string): Promise<Uint8Array>;
  runArchive(id: string): Promise<Uint8Array>;
  uploadInput(id: string, name: string, data: Blob): Promise<unknown>;
}

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
  datasets(): Promise<DatasetCatalog>;
  checkConfig(request: CheckConfigRequest): Promise<CheckConfigResult>;
  checkInputs(request: CheckInputsRequest): Promise<Parsed<typeof zInputFacts>>;
  listRuns(): Promise<Parsed<typeof zRunList>>;
  createRun(request: CreateRunRequest): Promise<Parsed<typeof zRunRecord>>;
  getRun(id: string): Promise<Parsed<typeof zRunRecord>>;
  startRun(id: string, request?: StartRunRequest): Promise<Parsed<typeof zRunRecord>>;
  updateRun(id: string, request: UpdateRunRequest): Promise<Parsed<typeof zRunSummary>>;
  cancelRun(id: string): Promise<boolean>;
  deleteRun(id: string): Promise<void>;
  restoreRun(id: string): Promise<Parsed<typeof zRunSummary>>;
  purgeRun(id: string): Promise<void>;
  followRun(id: string, options?: FollowRunOptions): Promise<TerminalEventResult>;
  runFiles(id: string): Promise<Array<Parsed<typeof zRunFile>>>;
  readRunFile(id: string, path: string): Promise<Uint8Array>;
  runArchive(id: string): Promise<Uint8Array>;
  uploadInput(id: string, name: string, data: Blob): Promise<Parsed<typeof zUploadedInput>>;
  timetree(config: TimetreeConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  optimize(config: OptimizeConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  prune(config: PruneConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  ancestral(config: AncestralConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  clock(config: ClockConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  mugration(config: MugrationConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
}

const zCancelResponse = z.object({ cancelled: z.boolean() });

export function createBridge(transport: BridgeTransport): TreeTimeBridge {
  async function followRun(id: string, options: FollowRunOptions = {}): Promise<TerminalEventResult> {
    let terminal: TerminalEventResult | undefined;

    const eventOptions: TransportEventOptions = {
      from: options.from ?? 0,
      onEvent: (data) => {
        const event = parseRunEvent(data);

        if (event.type === "terminal") {
          terminal = event.data;
        }

        options.onEvent?.(event);
      },
    };

    if (options.signal !== undefined) {
      eventOptions.signal = options.signal;
    }

    await transport.runEvents(id, eventOptions);

    if (terminal === undefined) {
      throw new RunEndedError(id);
    }

    return terminal;
  }

  async function run(
    command: AppCommand,
    config: unknown,
    options: CommandOptions = {},
  ): Promise<CommandOutcomeResult> {
    if (options.signal?.aborted === true) {
      throw new CancelledError();
    }

    const request: CreateRunRequest = { command, config, defer_start: false };

    if (options.title !== undefined) {
      request.title = options.title;
    }

    const record = zRunRecord.parse(await transport.createRun(request));

    const cancel = () => {
      void transport.cancelRun(record.id).catch((error: unknown) => {
        console.warn("[TreeTime] cancellation request failed", error);
      });
    };

    options.onStarted?.(record.id);
    options.signal?.addEventListener("abort", cancel);

    try {
      const terminal = await followRun(record.id, {
        onEvent: (event) => {
          dispatchEvent(event, options);
        },
      });

      return commandOutcome(terminal);
    } finally {
      options.signal?.removeEventListener("abort", cancel);
    }
  }

  return {
    async version() {
      return zVersionInfo.parse(await transport.version());
    },
    async datasets() {
      return zDatasetCatalog.parse(await transport.datasets());
    },
    async checkConfig(request) {
      return zCheckConfigResponse.parse(await transport.checkConfig(request));
    },
    async checkInputs(request) {
      return zInputFacts.parse(await transport.checkInputs(request));
    },
    async listRuns() {
      return zRunList.parse(await transport.listRuns());
    },
    async createRun(request) {
      return zRunRecord.parse(await transport.createRun(request));
    },
    async getRun(id) {
      return zRunRecord.parse(await transport.getRun(id));
    },
    async startRun(id, request = {}) {
      return zRunRecord.parse(await transport.startRun(id, request));
    },
    async updateRun(id, request) {
      return zRunSummary.parse(await transport.updateRun(id, request));
    },
    async cancelRun(id) {
      return zCancelResponse.parse(await transport.cancelRun(id)).cancelled;
    },
    deleteRun: (id) => transport.deleteRun(id),
    async restoreRun(id) {
      return zRunSummary.parse(await transport.restoreRun(id));
    },
    purgeRun: (id) => transport.purgeRun(id),
    followRun,
    async runFiles(id) {
      return z.array(zRunFile).parse(await transport.runFiles(id));
    },
    readRunFile: (id, path) => transport.readRunFile(id, path),
    runArchive: (id) => transport.runArchive(id),
    async uploadInput(id, name, data) {
      return zUploadedInput.parse(await transport.uploadInput(id, name, data));
    },
    timetree: (config, options) => run("timetree", config, options),
    optimize: (config, options) => run("optimize", config, options),
    prune: (config, options) => run("prune", config, options),
    ancestral: (config, options) => run("ancestral", config, options),
    clock: (config, options) => run("clock", config, options),
    mugration: (config, options) => run("mugration", config, options),
  };
}

function commandOutcome(terminal: TerminalEventResult): CommandOutcomeResult {
  if (terminal.status === "ok") {
    return terminal.result;
  }

  if (terminal.status === "cancelled") {
    throw new CancelledError();
  }

  if (terminal.status === "error") {
    throw new CommandError(terminal.job_id, terminal.message, terminal.causes);
  }

  throw new CommandError(terminal.job_id, "the run was interrupted because the process that ran it stopped", []);
}

function dispatchEvent(event: RunEventResult, options: CommandOptions): void {
  switch (event.type) {
    case "started":
      break;
    case "progress":
      options.onProgress?.(event.data);
      break;
    case "log":
      options.onLog?.(event.data);
      break;
    case "iteration":
      options.onIteration?.(event.data);
      break;
    case "terminal":
      break;
  }
}

export function parseRunEvent(data: unknown): RunEventResult {
  return zRunEvent.parse(data);
}
