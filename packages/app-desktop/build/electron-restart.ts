import { setTimeout as delay } from "node:timers/promises";

const EXIT_TIMEOUT_MS = 10_000;

export interface ElectronChild {
  readonly exitCode: number | null;
  readonly signalCode: NodeJS.Signals | null;
  removeAllListeners(event: "exit"): void;
  once(event: "exit", listener: () => void): void;
  kill(signal?: NodeJS.Signals): boolean;
}

export interface RestartOptions {
  current: () => ElectronChild | undefined;
  start: () => Promise<void>;
  exitDeadline?: () => Promise<void>;
}

export function restartChain(): (options: RestartOptions) => Promise<void> {
  let chain = Promise.resolve();

  return async (options) => {
    const next = chain.then(async () => restart(options));
    chain = next.catch(() => undefined);

    return next;
  };
}

async function restart({ current, start, exitDeadline = exitTimeout }: RestartOptions): Promise<void> {
  const child = current();

  if (child !== undefined && child.exitCode === null && child.signalCode === null) {
    await stop(child, exitDeadline);
  }

  await start();
}

async function stop(child: ElectronChild, exitDeadline: () => Promise<void>): Promise<void> {
  child.removeAllListeners("exit");

  const exited = new Promise<"exited">((resolve) => {
    child.once("exit", () => resolve("exited"));
  });

  child.kill();

  const deadline = exitDeadline().then(() => "deadline" as const);

  if ((await Promise.race([exited, deadline])) === "deadline") {
    child.kill("SIGKILL");
    await exited;
  }
}

async function exitTimeout(): Promise<void> {
  await delay(EXIT_TIMEOUT_MS);
}
