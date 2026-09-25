import * as z from "zod";

import type {
  AncestralConfig,
  AppCommand,
  AuspiceDocument,
  CheckConfigRequest,
  CladeRequest,
  CheckInputsRequest,
  ClockConfig,
  CreateRunRequest,
  DatasetCatalog,
  ErrorResponse,
  LogEvent,
  MugrationConfig,
  OptimizeConfig,
  ProgressEvent,
  PruneConfig,
  RunConfigRequest,
  StartRunRequest,
  TimetreeConfig,
  UpdateRunRequest,
  VersionInfo,
} from "./generated/types.gen";
import {
  zAuspiceDocument,
  zCancelRunResponse,
  zCheckConfigResponse,
  zCladeInRuns,
  zCommandOutcome,
  zDatasetCatalog,
  zDesktopRequest,
  zErrorResponse,
  zInputFacts,
  zIterationEvent,
  zRunConfigResponse,
  zRunEvent,
  zRunFile,
  zRunComparison,
  zRunList,
  zRunRecord,
  zRunResults,
  zRunSummary,
  zTerminalEvent,
  zUploadedInput,
  zVersionInfo,
} from "./generated/zod.gen";

export type Parsed<S extends z.ZodType> = z.infer<S>;

export type CheckConfigResult = Parsed<typeof zCheckConfigResponse>;

export type RunConfigResult = Parsed<typeof zRunConfigResponse>;

export type RunSummaryResult = Parsed<typeof zRunSummary>;

export type RunRecordResult = Parsed<typeof zRunRecord>;

export type InputFactsResult = Parsed<typeof zInputFacts>;

export type DesktopRequestInput = z.input<typeof zDesktopRequest>;

export type CheckConfigInput = Omit<CheckConfigRequest, "input_facts"> & { input_facts?: InputFactsResult | null };

export type RunResultsResult = Parsed<typeof zRunResults>;

export type RunComparisonResult = Parsed<typeof zRunComparison>;

export type CladeInRunsResult = Parsed<typeof zCladeInRuns>;

export type RunEventResult = Parsed<typeof zRunEvent>;

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
  onEvent: (event: unknown) => RunEventResult;
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

export class BridgeError extends Error {
  readonly response: ErrorResponse;

  constructor(response: ErrorResponse, context?: string) {
    const chain = [response.message, ...response.causes].join(": ");
    super(context === undefined ? chain : `${context}: ${chain}`);
    this.name = "BridgeError";
    this.response = response;
  }
}

export function bridgeErrorFromText(text: string, context?: string): BridgeError {
  const parsed = zErrorResponse.safeParse(parseJson(text));
  const response: ErrorResponse = parsed.success ? parsed.data : { code: "internal_error", message: text, causes: [] };

  return new BridgeError(response, context);
}

function parseJson(text: string): unknown {
  try {
    const value: unknown = JSON.parse(text);

    return value;
  } catch {
    return undefined;
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
  checkConfig(request: CheckConfigInput): Promise<unknown>;
  runConfig(request: RunConfigRequest): Promise<unknown>;
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
  saveRunFile(id: string, path: string, name: string): Promise<boolean>;
  saveRunArchive(id: string, name: string): Promise<boolean>;
  uploadInput(id: string, name: string, data: Blob): Promise<unknown>;
  runResults(id: string): Promise<unknown>;
  runAuspice(id: string): Promise<unknown>;
  compareRuns(id: string, other: string): Promise<unknown>;
  cladeInRuns(request: CladeRequest): Promise<unknown>;
}

export interface TreeTimeBridge {
  version(): Promise<VersionInfo>;
  datasets(): Promise<DatasetCatalog>;
  checkConfig(request: CheckConfigInput): Promise<CheckConfigResult>;
  runConfig(request: RunConfigRequest): Promise<RunConfigResult>;
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
  saveRunFile(id: string, path: string, name: string): Promise<boolean>;
  saveRunArchive(id: string, name: string): Promise<boolean>;
  uploadInput(id: string, name: string, data: Blob): Promise<Parsed<typeof zUploadedInput>>;
  runResults(id: string): Promise<RunResultsResult>;
  runAuspice(id: string): Promise<AuspiceDocument>;
  compareRuns(id: string, other: string): Promise<RunComparisonResult>;
  cladeInRuns(request: CladeRequest): Promise<CladeInRunsResult>;
  timetree(config: TimetreeConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  optimize(config: OptimizeConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  prune(config: PruneConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  ancestral(config: AncestralConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  clock(config: ClockConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
  mugration(config: MugrationConfig, options?: CommandOptions): Promise<CommandOutcomeResult>;
}

export function createBridge(transport: BridgeTransport): TreeTimeBridge {
  async function followRun(id: string, options: FollowRunOptions = {}): Promise<TerminalEventResult> {
    let terminal: TerminalEventResult | undefined;

    const eventOptions: TransportEventOptions = {
      from: options.from ?? 0,
      onEvent: (data) => {
        const event = parseRunEvent(data);

        if (event.type === "log") {
          logToConsole(event.data);
        }

        if (event.type === "terminal") {
          terminal = event.data;
        }

        options.onEvent?.(event);

        return event;
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
    async runConfig(request) {
      return zRunConfigResponse.parse(await transport.runConfig(request));
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
      return zCancelRunResponse.parse(await transport.cancelRun(id)).cancelled;
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
    saveRunFile: (id, path, name) => transport.saveRunFile(id, path, name),
    saveRunArchive: (id, name) => transport.saveRunArchive(id, name),
    async uploadInput(id, name, data) {
      return zUploadedInput.parse(await transport.uploadInput(id, name, data));
    },
    async runResults(id) {
      return zRunResults.parse(await transport.runResults(id));
    },
    async runAuspice(id) {
      return zAuspiceDocument.parse(await transport.runAuspice(id));
    },
    async compareRuns(id, other) {
      return zRunComparison.parse(await transport.compareRuns(id, other));
    },
    async cladeInRuns(request) {
      return zCladeInRuns.parse(await transport.cladeInRuns(request));
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

function logToConsole(log: LogEvent): void {
  switch (log.level) {
    case "error":
      console.error(`[TreeTime] ${log.message}`);
      break;
    case "warn":
      console.warn(`[TreeTime] ${log.message}`);
      break;
    case "info":
    case "debug":
    case "trace":
      console.log(`[TreeTime] [${log.level}] ${log.message}`);
      break;
  }
}

export function parseRunEvent(data: unknown): RunEventResult {
  return zRunEvent.parse(data);
}
