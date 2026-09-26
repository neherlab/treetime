import { CancelledError, CommandError, RunEndedError } from "@neherlab/app-contracts";
import type { CommandOutcome, CreateRunRequest, TerminalEvent } from "@neherlab/app-contracts";
import { runsCancel, runsCreate, type ApiClient } from "@neherlab/app-contracts/client";

import type { RunEvent } from "../results/progress";
import { runEventStream, type RunEventStreamOptions } from "./events";

export interface FollowRunOptions extends RunEventStreamOptions {
  onEvent?: ((event: RunEvent) => void) | undefined;
}

export type CommandRequest = Omit<CreateRunRequest, "defer_start">;

export interface RunCommandOptions extends Omit<FollowRunOptions, "from" | "signal"> {
  signal?: AbortSignal | undefined;
  onStarted?: ((id: string) => void) | undefined;
}

export async function followRun(client: ApiClient, id: string, options: FollowRunOptions = {}): Promise<TerminalEvent> {
  const { onEvent, ...streamOptions } = options;

  for await (const event of runEventStream(client, id, streamOptions)) {
    onEvent?.(event);

    if (event.type === "terminal") {
      return event.data;
    }
  }

  throw options.signal?.aborted === true ? new CancelledError() : new RunEndedError(id);
}

export async function runCommand(
  client: ApiClient,
  request: CommandRequest,
  options: RunCommandOptions = {},
): Promise<CommandOutcome> {
  const { signal, onStarted, ...followOptions } = options;

  if (signal?.aborted === true) {
    throw new CancelledError();
  }

  const { data: record } = await runsCreate({ client, body: { ...request, defer_start: false }, throwOnError: true });

  const cancel = () => {
    void cancelRun(client, record.id);
  };

  onStarted?.(record.id);
  signal?.addEventListener("abort", cancel, { once: true });

  try {
    return commandOutcome(await followRun(client, record.id, followOptions));
  } finally {
    signal?.removeEventListener("abort", cancel);
  }
}

async function cancelRun(client: ApiClient, id: string): Promise<void> {
  try {
    await runsCancel({ client, path: { id }, throwOnError: true });
  } catch (error: unknown) {
    console.warn("[TreeTime] cancellation request failed", error);
  }
}

function commandOutcome(terminal: TerminalEvent): CommandOutcome {
  if (terminal.status === "ok") {
    return terminal.result;
  }

  throw terminalError(terminal);
}

function terminalError(terminal: Exclude<TerminalEvent, { status: "ok" }>): Error {
  switch (terminal.status) {
    case "cancelled":
      return new CancelledError();
    case "error":
      return new CommandError(terminal.job_id, terminal.message, terminal.causes);
    case "interrupted":
      return new CommandError(terminal.job_id, "the run was interrupted because the process that ran it stopped", []);
    default:
      return new Error(`unknown terminal status of ${JSON.stringify(terminal satisfies never)}`);
  }
}
